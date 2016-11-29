__author__ = 'churas'

import os
import logging

from d3r.celpp.task import D3RTask
from d3r.celpp.task import D3RParameters
from d3r.celpp.proteinligprep import ProteinLigPrepTask
from d3r.celpp.dataimport import DataImportTask
from d3r.celpp.participant import ParticipantDatabaseFromCSVFactory
from d3r.celpp import util
from d3r.celpp.task import SmtpEmailer
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

    def _get_value_of_replytoaddress(self):
        """Attempts to get value of replytoaddress from `get_args()`
        :returns None if replytoaddress is None or does not exist in object
        """
        try:
            return self.get_args().replytoaddress
        except AttributeError:
            logger.debug('replytoaddress not set in get_args()')
            return None

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
        emailer = EvaluationEmailer(participant_db,
                                    self._get_value_of_replytoaddress())

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

        return scoring_tasks


class EvaluationEmailer(object):
    """Sends evaluation email
    """
    def __init__(self, participant_database, reply_to_address,
                 smtp='localhost', smtpport=25):
        self._participantdatabase = participant_database
        self._alt_smtp_emailer = None
        self._msg_log = None
        self._reply_to_address = reply_to_address
        self._smtp = smtp
        self._smtpport = smtpport

    def set_alternate_smtp_emailer(self, emailer):
        """Sets alternate smtp emailer
        """
        self._alt_smtp_emailer = emailer

    def _get_smtp_emailer(self):
        """Gets `SmtpEmailer` object or alternate
        """
        if self._alt_smtp_emailer is not None:
            return self._alt_smtp_emailer
        return SmtpEmailer(smtp_host=self._smtp,
                           port=self._smtpport)

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

    def _get_reply_to_address(self, from_address):
        """Gets reply to address
        :returns: reply to email address
        """
        if self._reply_to_address is None:
            return from_address
        return self._reply_to_address

    def _generate_external_submission_email_body(self, etask):
        """Creates body of email and subject
        :returns subject,body: as strings
        """
        weekno = util.get_celpp_week_number_from_path(etask.get_path())

        msg = 'Dear CELPP Participant,\n\nHere are your docking ' \
              'evaluation results (RMSD, Angstroms) for CELPP week ' +\
              str(weekno) + '\n\n'
        msg += etask.get_evaluation_summary()

        msg += '\n\nSincerely,\n\nCELPP Automation'

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
            emailer = self._get_smtp_emailer()

            to_list = self._get_external_submitter_email(etask)
            if to_list is None:
                logger.debug('No external submitter email, just returning')
                return

            subject, msg = self\
                ._generate_external_submission_email_body(etask)

            from_addr = emailer.generate_from_address_using_login_and_host()

            reply_to = self._get_reply_to_address(from_addr)

            rmsd = Attachment(etask.get_rmsd_txt(), 'rmsd.txt')

            emailer.send_email(from_addr, to_list, subject, msg,
                               reply_to=reply_to,
                               attachments=[rmsd])
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
    RMSD_PICKLE = 'RMSD.pickle'
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

    def get_rmsd_pickle(self):
        """Returns full path to RMSD.pickle file
        :returns: full path to RMSD.pickle file
        """
        return os.path.join(self.get_dir(), EvaluationTask.RMSD_PICKLE)

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
            logger.error('Parent DockTask is None, so guid is'
                         ' NOT extractable')
            return None

        return self._docktask.get_name()\
            .replace(EvaluationTask.EXT_SUBMISSION_SUFFIX, '')

    def get_uploadable_files(self):
        """Returns list of files that can be uploaded to remote server

           List will contain these files if found on the file system
           plus stderr/stdout files

           RMSD.txt
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

            rmsd = self.get_rmsd_txt()
            if os.path.isfile(rmsd):
                file_list.append(rmsd)

            rmsdpickle = self.get_rmsd_pickle()
            if os.path.isfile(rmsdpickle):
                file_list.append(rmsdpickle)

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

        #
        # --pdbdb <path to pdb.extracted> --dockdir <stage.4.glide> \
        # --outdir <path to stage.5.glide.evaluation>
        #
        cmd_to_run = (self.get_args().evaluation + ' --pdbdb ' +
                      self.get_args().pdbdb + ' --dockdir ' +
                      self._docktask.get_dir() +
                      ' --outdir ' + self.get_dir())

        eval_name = os.path.basename(self.get_args().evaluation)

        self.run_external_command(eval_name, cmd_to_run,
                                  True)

        # attempt to send evaluation email
        try:
            self._emailer.send_evaluation_email(self)
            self.append_to_email_log(self._emailer.get_message_log())
        except Exception as e:
            logger.exception('Caught exception trying to send evaluation '
                             'email')
            self.append_to_email_log('Caught exception trying to send '
                                     'evaluation email ' + str(e) + '\n')

        # assess the result
        self.end()
