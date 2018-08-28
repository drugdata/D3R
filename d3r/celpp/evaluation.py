__author__ = 'churas'

import os
import logging
import tarfile
import re
import json
import pickle
import requests
from requests.auth import HTTPBasicAuth

from d3r.celpp.task import SmtpEmailerFactory
from d3r.celpp.task import D3RTask
from d3r.celpp.task import WebsiteServiceConfig
from d3r.celpp.task import D3RParameters
from d3r.celpp.proteinligprep import ProteinLigPrepTask
from d3r.celpp.dataimport import DataImportTask
from d3r.celpp.blastnfilter import BlastNFilterTask
from d3r.celpp.challengedata import ChallengeDataTask
from d3r.celpp.participant import ParticipantDatabaseFromCSVFactory
from d3r.celpp import util
from d3r.celpp.task import Attachment


logger = logging.getLogger(__name__)


class PathNotDirectoryError(Exception):
    """Path is not a directory Error
    """
    pass


class EvaluationTaskFactory(object):
    """Factory class to generate ScoringTask objects

       This factory examines a celpp week directory for
       all docking tasks.  The code then generates
       ScoringTask objects for all eligible docking tasks
    """
    prep = ProteinLigPrepTask('/foo', D3RParameters())
    DOCKSTAGE = prep.get_stage() + 1
    DOCKSTAGE_PREFIX = (D3RTask.STAGE_DIRNAME_PREFIX + '.' +
                        str(DOCKSTAGE) + '.')
    SCORING_SUFFIX = 'evaluation'
    WEB_DATA_SUFFIX = 'webdata'

    def __init__(self, path, theargs):
        """Constructor
        """
        self.set_path(path)
        self.set_args(theargs)

    def set_args(self, theargs):
        """ Sets args
        :param theargs: arguments to set
        """
        self._args = theargs

    def get_args(self):
        """Gets args passed into constructor or via set_args()
        """
        return self._args

    def set_path(self, path):
        """Sets path used to look for docking tasks
        """
        self._path = path

    def get_path(self):
        """Gets path used to look for docking tasks
        """
        return self._path

    def _get_participant_database(self):
        """Creates `ParticipantDatabase`
        :returns: ParticipantDatabase
        """
        dimport = DataImportTask(self.get_path(), self.get_args())
        csvfile = dimport.get_participant_list_csv()
        pfac = ParticipantDatabaseFromCSVFactory(csvfile)
        return pfac.get_participant_database()

    def _update_priorities_of_tasks(self, etasks,
                                    participant_db):
        """Updates priorities for EvaluationTask objects in
           list of etask objects from values obtained in
           participant_db ParticipantDatabase
           :param etasks: list of EvaluationTask objects, it is assumed to
                          never be None
           :param participant_db: ParticipantDatabase object, assumed never
                                  to be none
           :returns: same list of EvaluationTask objects with priority set
        """
        for task in etasks:
            guid = task.get_guid_for_task()
            if guid is not None:
                logger.debug('Looking for participant with guid: ' + guid)
                p = participant_db.get_participant_by_guid(guid)
                if p is not None:
                    task.set_priority(p.get_priority())
                    logger.debug('Setting priority for ' + guid +
                                 ' to ' + str(task.get_priority()))
                else:
                    logger.debug('No participant found to match guid: ' + guid)
        return etasks

    def _sort_tasks_by_participant_priority(self, etasks,
                                            participant_db):
        """Sorts `EvaluationTask` objects in `etasks` by priority
           set for participants in `participant_db` The sorting
           goes as follows. Participants with highest get_priority()
           value go first, identical priority ordering is arbitrary
           and any EvaluationTasks without priority are put at end
           of list in arbitrary order
        :returns: list of sorted EvaluationTask objects
        """
        if etasks is None:
            logger.debug('No EvaluationTasks to sort')
            return etasks

        if participant_db is None:
            logger.warning('Participant Database is None, cannot sort')
            return etasks

        # update priorities of etasks
        updatedtasks = self._update_priorities_of_tasks(etasks, participant_db)
        updatedtasks.sort(reverse=True, key=lambda task: task.get_priority())
        return updatedtasks

    def get_evaluation_tasks(self):
        """Generate EvaluationTasks

           This method examines the path directory
           set via the constructor or set_path() method
           for all stage 4 tasks excluding tasks
           that end with 'webdata'  A EvaluationTask
           object is created for each of these tasks
           and returned in a list.
           :return: list of EvaluationTask objects or empty list if none found
        """
        path = self.get_path()
        logger.debug('Examining ' + path + ' for docking tasks')
        if not os.path.isdir(path):
            raise PathNotDirectoryError(path + ' is not a directory')
        scoring_tasks = []

        path_list = os.listdir(path)

        participant_db = self._get_participant_database()
        efac = SmtpEmailerFactory(self.get_args())
        emailer = EvaluationEmailer(participant_db,
                                    efac.get_smtp_emailer())

        for entry in path_list:
            logger.debug('Checking if ' + entry + ' is a docking task')
            full_path = os.path.join(path, entry)
            if os.path.isdir(full_path):
                if entry.startswith(EvaluationTaskFactory.DOCKSTAGE_PREFIX):
                    if entry.endswith(EvaluationTaskFactory.WEB_DATA_SUFFIX):
                        logger.debug('Skipping ' + entry + ' due to suffix')
                        continue

                    # we have a valid docking path
                    docktask = D3RTask(path, self.get_args())
                    docktask.set_stage(EvaluationTaskFactory.DOCKSTAGE)
                    docktask.set_name(entry[
                                      len(EvaluationTaskFactory
                                          .DOCKSTAGE_PREFIX):])
                    stask = EvaluationTask(path,
                                           docktask.get_name() + '.' +
                                           EvaluationTaskFactory.
                                           SCORING_SUFFIX,
                                           docktask,
                                           self.get_args())
                    if stask.can_run():
                        logger.debug('Adding task ' + stask.get_name())
                        stask.set_evaluation_emailer(emailer)
                        scoring_tasks.append(stask)
                    else:
                        if stask.get_error() is None:
                            logger.debug(stask.get_name() + ' cannot be' +
                                         ' added, no error though')
                        else:
                            logger.debug(stask.get_name() + ' cannot be' +
                                         ' added : ' + stask.get_error())

        return self._sort_tasks_by_participant_priority(scoring_tasks,
                                                        participant_db)


class EvaluationEmailer(object):
    """Sends evaluation email
    """
    def __init__(self, participant_database, emailer):
        self._participantdatabase = participant_database
        self._emailer = emailer
        self._msg_log = None

    def _append_to_message_log(self, msg):
        if self._msg_log is None:
            self._msg_log = str(msg)
            return
        self._msg_log += str(msg)

    def get_message_log(self):
        """Gets the error log for this object
        :returns: Error log as string
        """
        return self._msg_log

    def _get_external_submitter_email(self, etask):
        """Extracts external submitter email from ParticipantDatabase

        :return: list of external submitter email or emails ot None if
        not found
        """
        # need to get guid from task name stage.#.GUID.extsubmission
        guid = etask.get_guid_for_task()
        if guid is None:
            logger.error('guid is None')
            self._append_to_message_log('\nUnable to extract guid\n')
            return None

        # then need to get email for guid
        part = self._participantdatabase.get_participant_by_guid(guid)

        if part is None:
            logger.error('No participant found with guid: ' + guid)
            self._append_to_message_log('\nNo participant found with guid: ' +
                                        guid + '\n')
            return None

        if part.get_email() is None:
            logger.error('Email not set for participant')
            self._append_to_message_log('\nEmail not set for participant\n')
            return None

        return part.get_email().split(',')

    def _generate_external_submission_email_body(self, etask):
        """Creates body of email and subject
        :returns subject,body: as strings
        """
        weekno = util.get_celpp_week_number_from_path(etask.get_path())

        msg = 'Dear CELPP Participant,\n\nHere are your docking ' \
              'evaluation results as RMSD in Angstroms for CELPP week ' +\
              str(weekno) + '\n\n'
        msg += ('Note: The value in (parentheses) by each RMSD is the ' +
                'distance, in Angstroms, between submitted ligand center ' +
                'of mass and the crystal ligand center of mass. The final '
                'column is the distance between CELPP-provided pocket '
                'center and the crystal ligand center of mass.\n\n')
        msg += etask.get_evaluation_summary()

        msg += ('\n\nSincerely,\n\nCELPP Automation ' +
                etask.get_program_version())

        guid = etask.get_guid_for_task()
        subject = (D3RTask.SUBJECT_LINE_PREFIX +
                   'Week ' + str(weekno) + ' evaluation results for ' +
                   str(guid))
        return subject, msg

    def send_evaluation_email(self, etask):
        """Sends evaluation email appending issues to message log
        """

        # clear message log fix for issue #99
        self._msg_log = None
        if etask is None:
            logger.error('Task passed in is None')
            self._append_to_message_log('\nTask passed in is None\n')
            return

        if etask.is_external_submission() is False:
            logger.debug('Not an external submission, just returning')
            self._append_to_message_log('\nNot an external submission\n')
            return

        if self._participantdatabase is None:
            logger.error('Participant Database is None')
            self._append_to_message_log('\nParticipant database is None '
                                        'cannot send'
                                        ' docking evaluation email!!!\n')
            return
        try:

            to_list = self._get_external_submitter_email(etask)
            if to_list is None:
                logger.debug('No external submitter email, just returning')
                return

            subject, msg = self\
                ._generate_external_submission_email_body(etask)

            rmsd = Attachment(etask.get_rmsd_txt(), 'rmsd.txt')

            self._emailer.send_email(to_list, subject, msg, attachments=[rmsd])
            self._append_to_message_log('\nSent evaluation email to: ' +
                                        ", ".join(to_list) + '\n')
        except Exception as e:
            logger.exception('Caught exception')
            self._append_to_message_log('\nCaught exception trying to email '
                                        'participant : ' + str(e) + '\n')


class EvaluationTask(D3RTask):
    """Performs evaluation

    """

    FINAL_LOG = 'final.log'
    RMSD_TXT = 'RMSD.txt'
    RMSD_JSON = 'RMSD.json'
    RMSD_PICKLE = 'RMSD.pickle'
    RMSD_CSV = 'RMSD.csv'

    EVAL_EXITCODEFILE = 'evaluate.exitcode'

    EXT_SUBMISSION_SUFFIX = '.extsubmission'
    SCORE_DIR = 'score'
    COMPLEX_SUFFIX = '_complex.pdb'

    PDB_FILES = [SCORE_DIR + os.sep + 'crystal.pdb']

    def __init__(self, path, name, docktask, args):
        super(EvaluationTask, self).__init__(path, args)
        self.set_name(name)
        self.set_stage(EvaluationTaskFactory.DOCKSTAGE + 1)
        self.set_status(D3RTask.UNKNOWN_STATUS)
        self._docktask = docktask
        self._emailer = None
        self._priority = 0
        self._webserviceconfig = None

        try:
            config = self.get_args().websiteserviceconfig
            self._webserviceconfig = WebsiteServiceConfig(configfile=config)
        except AttributeError:
            pass
        try:
            self._week_num = util.\
                get_celpp_week_number_from_path(self.get_path())
            self._year = util.get_celpp_year_from_path(self.get_dir())
        except Exception:
            self._week_num = 0
            self._year = 0

    def set_priority(self, priority):
        """Sets priority of when evaluation task should be run
        :param priority: int denoting priority, higher values mean
                         higher priority and 0 is usually default
        """
        self._priority = priority

    def get_priority(self):
        """Gets priority of this EvaluationTask
        :returns: int denoting priority
        """
        return self._priority

    def _write_evaluate_exitcode_file(self, ecode):
        """Writes evaluate.py exit code to
           self.get_dir().EVAL_EXITCODEFILE
           If there is an error of any kind this method will
           catch the exception and log it to email and to
           the logger
        """
        try:
            ecodefile = os.path.join(self.get_dir(),
                                     EvaluationTask.EVAL_EXITCODEFILE)
            logger.debug('Attempting to write to ' + str(ecodefile))
            with open(ecodefile, 'w') as f:
                f.write(str(ecode) + '\n')
                f.flush()
        except Exception as e:
            logger.exception('Error writing exit code file')
            self.append_to_email_log('Caught exception trying to '
                                     'write exit code file ' + str(e) + '\n')

    def set_evaluation_emailer(self, emailer):
        self._emailer = emailer

    def is_external_submission(self):
        """Checks if this evaluation task is analyzing an external submission
        The check is done by looking at suffix of name to see if it matches
        `ExternalDataSubmissionTask.EXT_SUBMISSION_SUFFIX`
        :returns: True for yes otherwise False
        """
        if self._docktask is None:
            return False

        if self._docktask.get_name().endswith(
                EvaluationTask.EXT_SUBMISSION_SUFFIX):
            return True
        return False

    def get_rmsd_txt(self):
        """Returns full path to RMSD.txt file
        :returns: full path to RMSD.txt file
        """
        return os.path.join(self.get_dir(), EvaluationTask.RMSD_TXT)

    def get_rmsd_json(self):
        """Returns full path to RMSD.json file
        :returns: full path to RMSD.json file
        """
        return os.path.join(self.get_dir(), EvaluationTask.RMSD_JSON)

    def get_rmsd_pickle(self):
        """Returns full path to RMSD.pickle file
        :returns: full path to RMSD.pickle file
        """
        return os.path.join(self.get_dir(), EvaluationTask.RMSD_PICKLE)

    def get_rmsd_csv(self):
        """Returns full path to RMSD.csv file
        :returns: full path to RMSD.csv file
        """
        return os.path.join(self.get_dir(), EvaluationTask.RMSD_CSV)

    def get_evaluation_summary(self):
        """Parses RMSD.txt and generates human readable summary
        evaluating docking
        :returns: report of docking as string.
        """
        rmsd = self.get_rmsd_txt()

        start_line = '\nEvaluation of docking\n=====================\n'

        if not os.path.isfile(rmsd):
            return start_line + 'No ' + rmsd + ' file found.\n'

        try:
            f = open(rmsd, 'r')
            summary = f.read()
            f.close()
            return start_line + summary + '\n'
        except Exception as e:
            logger.exception('Caught exception')
            return (start_line + 'Unable to generate evaluation summary (' +
                    str(e) + ')\n')

    def get_guid_for_task(self):
        """Examines task name removing external submission suffix
           to hopefully obtain guid
        :returns: guid from task name
        """
        if self._docktask is None:
            logger.debug('Dock task is None going to strip off .* suffix')
            if self.get_name() is None:
                logger.error('Task name is none, giving up trying to get guid')
                return None
            return re.sub('\..*$', '', self.get_name())

        return self._docktask.get_name()\
            .replace(EvaluationTask.EXT_SUBMISSION_SUFFIX, '')

    def get_uploadable_files(self):
        """Returns list of files that can be uploaded to remote server

           List will contain these files if found on the file system
           plus stderr/stdout files

           RMSD.txt
           RMSD.json
           RMSD.pickle
           final.log
           pbdid/score/crystal.pdb
           pbdid/score/*_complex.pdb

           :returns: list of files that can be uploaded
        """
        # get the stderr/stdout files
        file_list = super(EvaluationTask, self).get_uploadable_files()

        out_dir = self.get_dir()

        try:
            final_log = os.path.join(out_dir, EvaluationTask.FINAL_LOG)
            if os.path.isfile(final_log):
                file_list.append(final_log)

            efile = os.path.join(out_dir, EvaluationTask.EVAL_EXITCODEFILE)
            if os.path.isfile(efile):
                file_list.append(efile)

            rmsd = self.get_rmsd_txt()
            if os.path.isfile(rmsd):
                file_list.append(rmsd)

            rmsdjson = self.get_rmsd_json()
            if os.path.isfile(rmsdjson):
                file_list.append(rmsdjson)

            rmsdpickle = self.get_rmsd_pickle()
            if os.path.isfile(rmsdpickle):
                file_list.append(rmsdpickle)

            rmsdcsv = self.get_rmsd_csv()
            if os.path.isfile(rmsdcsv):
                file_list.append(rmsdcsv)

            for entry in os.listdir(out_dir):
                full_path = os.path.join(out_dir, entry)
                if not os.path.isdir(full_path):
                    continue

                logger.debug('Looking for pdb files in ' + full_path)
                for pdb_name in EvaluationTask.PDB_FILES:
                    pdb = os.path.join(full_path, pdb_name)
                    if os.path.isfile(pdb):
                        file_list.append(pdb)
                score_dir = os.path.join(full_path,
                                         EvaluationTask.SCORE_DIR)
                if os.path.isdir(score_dir):
                    for score_entry in os.listdir(score_dir):
                        if score_entry.endswith(EvaluationTask.COMPLEX_SUFFIX):
                            pdb_f_path = os.path.join(score_dir, score_entry)
                            if os.path.isfile(pdb_f_path):
                                file_list.append(pdb_f_path)
        except OSError:
            logger.exception('Caught exception looking for pbdid folders')
        return file_list

    def can_run(self):
        """Determines if task can actually run

           This method first verifies the docking task
           has `D3RTask.COMPLETE_STATUS` for
           status.  The method then verifies this task does
           not already exist.  If above is not true then self.set_error()
           is set with information about the issue
           :return: True if can run otherwise False
        """
        self._can_run = False
        self._error = None
        self._docktask.update_status_from_filesystem()
        if self._docktask.get_status() != D3RTask.COMPLETE_STATUS:
            logger.info('Cannot run ' + self.get_name() + ' task ' +
                        'because ' + self._docktask.get_name() + ' task' +
                        ' has a status of ' + self._docktask.get_status())
            self.set_error(self._docktask.get_name() + ' task has ' +
                           self._docktask.get_status() + ' status')
            return False

        # check blast
        blastnfilter = BlastNFilterTask(self._path, self._args)
        blastnfilter.update_status_from_filesystem()
        if blastnfilter.get_status() != D3RTask.COMPLETE_STATUS:
            logger.info('Cannot run ' + self.get_name() + 'task ' +
                        'because ' + blastnfilter.get_name() + 'task' +
                        'has a status of ' + blastnfilter.get_status())
            self.set_error(blastnfilter.get_name() + ' task has ' +
                           blastnfilter.get_status() + ' status')
            return False

        # check this task is not complete and does not exist

        self.update_status_from_filesystem()
        if self.get_status() == D3RTask.COMPLETE_STATUS:
            logger.debug("No work needed for " + self.get_name() +
                         " task")
            return False

        if self.get_status() != D3RTask.NOTFOUND_STATUS:
            logger.warning(self.get_name() + " task was already " +
                           "attempted, but there was a problem")
            self.set_error(self.get_dir_name() + ' already exists and ' +
                           'status is ' + self.get_status())
            return False
        self._can_run = True
        return True

    def get_evaluationresult_filename(self, nosuffix=True):
        """Gets the evaluationresult filename
        :param nosuffix: boolean if set to True .tar.gz suffix is omitted
        :returns: Just filename as string in format
                  celpp_week<week no>_<year>_evalresults_<guid>
        """
        res = ('celpp_week' + str(self._week_num) + '_' + str(self._year) +
               '_evalresults_' + str(self.get_guid_for_task()))
        if nosuffix is False:
            res += '.tar.gz'
        return res

    def _create_evaluationresult_tarfile(self, evalresultdirname):
        """Creates a gzipped tarfile of the results of the evaluation
           that can be uploaded to remote server for download by
           participant. Tarfile will have name:
           `evalresultdirname`.tar.gz

           with contents under directory with name of file minus .tar.gz:

           RMSD.csv
           RMSD.txt
           RMSD.json
           start
           complete
           error
           evaluate.exitcode
           final.log

           :returns: string that is full path to tarfile.
        """
        evalfile = os.path.join(self.get_dir(), evalresultdirname + '.tar.gz')

        tar = tarfile.open(evalfile, 'w:gz')

        file_list = [EvaluationTask.FINAL_LOG,
                     EvaluationTask.EVAL_EXITCODEFILE,
                     EvaluationTask.RMSD_CSV, EvaluationTask.RMSD_JSON,
                     EvaluationTask.RMSD_JSON, EvaluationTask.COMPLETE_FILE,
                     EvaluationTask.ERROR_FILE, EvaluationTask.START_FILE]

        for entry in file_list:
            fullpath_entry = os.path.join(self.get_dir(), entry)
            if os.path.isfile(fullpath_entry):
                tar.add(fullpath_entry, arcname=evalresultdirname + '/' +
                        entry)
        tar.close()
        return evalfile

    def _upload_evaluationresult_file(self, evalresult_file):
        """Uploads `evalresult_file` to remote server
        """
        if evalresult_file is None:
            logger.error('evalresult_file is None')
            self.append_to_email_log('Skipping upload of evaluationresult '
                                     'since its None\n')
            return

        if not os.path.isfile(evalresult_file):
            logger.error(evalresult_file + ' file not found')
            self.append_to_email_log(evalresult_file + ' file not found '
                                                       'skipping upload\n')
            return

        uploader = self.get_file_transfer()
        if uploader is None:
            logger.warning('No uploader available to upload evaluation result '
                           'file')
            self.append_to_email_log('No uploader available to upload '
                                     'evaluation result file\n')
            return
        try:
            remote_dir = uploader.get_remote_evaluationresult_dir()
            if remote_dir is None:
                logger.warning('No remote evaluation result directory set for '
                               'ftp upload')
                self.append_to_email_log('No remote evaluation result '
                                         'directory set ' +
                                         'for ftp upload\n')
                return

            guid = self.get_guid_for_task()
            if guid is None:
                logger.error('Unable to get guid for task')
                self.append_to_email_log('Unable to get guid for task '
                                         'for ftp upload\n')
                return

            remote_dir = os.path.join(remote_dir, guid)
            logger.debug('Attempting to upload ' + evalresult_file + ' to '
                         + remote_dir)

            logger.debug('Connecting to remote server to upload challenge data'
                         'package')
            # connect to remote server
            uploader.connect()

            if uploader.upload_file_direct(evalresult_file, remote_dir,
                                           os.path.basename(evalresult_file)) \
                    is False:
                raise Exception(uploader.get_error_msg())

            self.append_to_email_log('\nEvaluation result tarball uploaded: \n'
                                     + uploader.get_upload_summary() + '\n')
        finally:
            uploader.disconnect()

    def _get_rmsd(self):
        """Loads RMSD for this task by parsing first the RMSD.json file
           and falling back to RMSD.pickle
           :returns output of json.load(RMSD.json) or pickle.load(RMSD.pickle)
        """
        rmsdjson = self.get_rmsd_json()
        rmsdpickle = self.get_rmsd_pickle()
        if os.path.isfile(rmsdjson):
            with open(rmsdjson, 'r') as fp:
                return json.load(fp)
        elif os.path.isfile(rmsdpickle):
            with open(rmsdpickle, 'r') as fp:
                return pickle.load(fp)
        else:
            return None

    def generate_rmsd_object(self):
        """Generates a dictionary object parsed from either RMSD.json or as
           a fallback as RMSD.pickle file that contains scores from this
           evaluation. Format matches the one defined for REST service here:
           https://github.com/drugdata/D3R/wiki/D3R-Website-REST-service

           """
        rmsdobj = {}
        rmsd_val = self._get_rmsd()

        if rmsd_val is None:
            logger.debug('RMSD object is None returning')
            return None

        rmsdobj[WebsiteServiceConfig.JSON_KEY] = rmsd_val

        version = self.get_version_from_start_file()

        wsource = 'notset'
        pname = 'notset'
        if self._webserviceconfig is not None:
            wsource = self._webserviceconfig.get_source()
            pname = self._webserviceconfig.get_portal_name()

        rmsdobj['version'] = version
        rmsdobj['source'] = wsource
        rmsdobj['week'] = int(self._week_num)
        rmsdobj['year'] = int(self._year)
        rmsdobj['portal_name'] = pname
        rmsdobj['submission_folder'] = self.get_guid_for_task()
        logger.debug('Successfully generated RMSD object for website service')
        return rmsdobj

    def post_rmsd_to_websiteservice(self, rmsdobj):
        """Posts RMSD scores to evaluation service
        """

        if rmsdobj is None:
            logger.debug('RMSD object is None skipping post to websiteservice')
            return

        if self._webserviceconfig is None:
            logger.warning('No service information available to post '
                           'evaluation results')
            self.append_to_email_log('No website service configuration found '
                                     'to post evaluation results\n')
            return

        theheader = dict()

        theheader[WebsiteServiceConfig.API_KEY_KEY] = self.\
            _webserviceconfig.get_apikey()
        theheader[WebsiteServiceConfig.
                  CONTENT_TYPE_KEY] = WebsiteServiceConfig.CONTENT_TYPE_VAL
        bauth = None
        if self._webserviceconfig.get_basicauth_user() is not None:
            bauth = HTTPBasicAuth(self._webserviceconfig.get_basicauth_user(),
                                  self._webserviceconfig.
                                  get_basicauth_password())
        logger.debug('Posting RMSD object to ' +
                     str(self._webserviceconfig.get_rmsd_url()))
        r = requests.post(self._webserviceconfig.get_rmsd_url(),
                          headers=theheader, json=rmsdobj, auth=bauth,
                          timeout=self._webserviceconfig.get_timeout())

        if r.status_code != 200:
            msg = ('website REST service returned code: ' +
                   str(r.status_code) + ' with text: ' + str(r.text) +
                   ' and json: ' + str(r.json))
            logger.error(msg)
            self.append_to_email_log('\n' + msg + '\n')
        else:
            logger.debug('Post returned code 200. Success')
        return

    def run(self):
        """Runs EvaluationTask after verifying dock was good

           Method requires can_run() to be called before hand with
           successful outcome
           Otherwise method invokes D3RTask.start then this method
           creates a directory and invokes evaluation script.  Upon
           completion results are analyzed and success or error status
           is set appropriately and D3RTask.end is invoked
           """
        super(EvaluationTask, self).run()

        if self._can_run is False:
            logger.debug(
                self.get_dir_name() + ' cannot run cause _can_run flag '
                                      'is False')
            return

        try:
            logger.debug('evaluation set to ' +
                         self.get_args().evaluation)
        except AttributeError:
            self.set_error('evaluation not set')
            self.end()
            return

        try:
            logger.debug('pdbdb set to ' +
                         self.get_args().pdbdb)
        except AttributeError:
            self.set_error('pdbdb not set')
            self.end()
            return

        try:
            evaltimeout = self.get_args().evaluationtimeout
            logger.debug('Setting evaluation timeout to ' +
                         str(evaltimeout))
        except AttributeError:
            evaltimeout = None

        try:
            killdelay = self.get_args().evaluationtimeoutkilldelay
            logger.debug('Setting evaluation kill delay to ' +
                         str(killdelay))
        except AttributeError:
            killdelay = 60

        blastnfilter = BlastNFilterTask(self._path, self._args)
        challenge = ChallengeDataTask(self._path, self._args)
        challdir = os.path.join(challenge.get_dir(),
                                challenge.get_celpp_challenge_data_dir_name())
        #
        # --pdbdb <path to pdb.extracted> --dockdir <stage.4.glide> \
        # --outdir <path to stage.5.glide.evaluation>
        #
        cmd_to_run = (self.get_args().evaluation + ' --pdbdb ' +
                      self.get_args().pdbdb + ' --dockdir ' +
                      self._docktask.get_dir() +
                      ' --blastnfilterdir ' +
                      blastnfilter.get_dir() +
                      ' --challengedir ' + challdir +
                      ' --outdir ' + self.get_dir())

        eval_name = os.path.basename(self.get_args().evaluation)

        ecode = self.run_external_command(eval_name, cmd_to_run,
                                          False,
                                          timeout=evaltimeout,
                                          kill_delay=killdelay)

        # write out evaluate exit code file
        self._write_evaluate_exitcode_file(ecode)

        # attempt to send evaluation email
        try:
            self._emailer.send_evaluation_email(self)
            self.append_to_email_log(self._emailer.get_message_log())
        except Exception as e:
            logger.exception('Caught exception trying to send evaluation '
                             'email')
            self.append_to_email_log('Caught exception trying to send '
                                     'evaluation email ' + str(e) + '\n')

        # attempt to post evaluation results to website REST service
        try:
            self.post_rmsd_to_websiteservice(self.generate_rmsd_object())
        except Exception as e:
            logger.exception('Not a show stopper, but caught exception '
                             'trying to post results to website rest '
                             'service')
            self.append_to_email_log('\nCaught exception trying to post '
                                     'evaluation results to website rest '
                                     'service ' + str(e) + '\n')
        # assess the result
        self.end()
