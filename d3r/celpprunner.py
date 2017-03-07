#! /usr/bin/env python

import sys
import os
import argparse
import psutil
import logging
from datetime import date

import d3r
from d3r.celpp import util
from d3r.celpp.task import D3RParameters


# create logger
logger = logging.getLogger('d3r.celpprunner')
DEFAULT_LOG_LEVEL = 'ERROR'
p = D3RParameters()
p.loglevel = DEFAULT_LOG_LEVEL
util.setup_logging(p)


from d3r.celpp.blastnfilter import BlastNFilterTask
from d3r.celpp.proteinligprep import ProteinLigPrepTask
from d3r.celpp.dataimport import DataImportTask
from d3r.celpp.glide import GlideTask
from d3r.celpp.evaluation import EvaluationTaskFactory
from d3r.celpp.makeblastdb import MakeBlastDBTask
from d3r.celpp.vina import AutoDockVinaTask

from d3r.celpp.challengedata import ChallengeDataTask
from d3r.celpp.chimeraprep import ChimeraProteinLigPrepTask
from d3r.celpp.filetransfer import FtpFileTransfer
from d3r.celpp.extsubmission import ExternalDataSubmissionFactory

from lockfile.pidlockfile import PIDLockFile


CREATE_CHALLENGE = 'createchallenge'
CHIMERA_PREP = 'chimeraprep'


def _get_lock(theargs, stage):
    """Create lock file to prevent this process from running on same data.

       This uses ``PIDLockFile`` to create a pid lock file in latest_weekly
       directory named celprunner.<stage>.lockpid
       If pid exists it is assumed the lock is held otherwise lock
       is broken and recreated

       :param theargs: return value from argparse and should contain
                       theargs.latest_weekly should be set to path
       :param stage: set to stage that is being run
       :return: ``PIDLockFile`` upon success
       :raises: LockException: If there was a problem locking
       :raises: Exception: If valid pid lock file already exists
       """
    mylockfile = os.path.join(theargs.latest_weekly, "celpprunner." +
                              stage + ".lockpid")
    logger.debug("Looking for lock file: " + mylockfile)
    lock = PIDLockFile(mylockfile, timeout=10)

    if lock.i_am_locking():
        logger.debug("My process id" + str(lock.read_pid()) +
                     " had the lock so I am breaking")
        lock.break_lock()
        lock.acquire(timeout=10)
        return lock

    if lock.is_locked():
        logger.debug("Lock file exists checking pid")
        if psutil.pid_exists(lock.read_pid()):
            raise Exception("celpprunner with pid " +
                            str(lock.read_pid()) +
                            " is running")

    lock.break_lock()
    logger.info("Acquiring lock")
    lock.acquire(timeout=10)
    return lock


def set_andor_create_latest_weekly_parameter(theargs):
    """Looks at theargs parameters to get celpp week directory

       What this method does varies by values in theargs
       parameter.

       If theargs.createweekdir is set then code will determine
       current celp week ie (dataset.week.#) and create it under
       the theargs.celppdir.

       Otherwise code will find the latest
       celpp week dir under the theargs.celppdir directory and set
       theargs.latest_weekly to this path.


    """
    try:
        if theargs.createweekdir:
            celp_week = util.get_celpp_week_of_year_from_date(date.today())
            logger.debug('Request to create new directory ' +
                         os.path.join(theargs.celppdir, str(celp_week[1]),
                                      'dataset.week.'+str(celp_week[0])))
            util.create_celpp_week_dir(celp_week, theargs.celppdir)
    except AttributeError:
        pass

    try:
        if theargs.customweekdir:
            theargs.latest_weekly = theargs.celppdir
        else:
            theargs.latest_weekly = \
                util.find_latest_weekly_dataset(theargs.celppdir)
    except AttributeError:
        theargs.latest_weekly = \
            util.find_latest_weekly_dataset(theargs.celppdir)
    return theargs


def run_stages(theargs):
    """Runs all the stages set in theargs.stage parameter


       Examines theargs.stage and splits it by comma to get
       list of stages to run.  For each stage found a lock file
       is created and run_stage is invoked with theargs.latest_weekly set to
       the output of util.find_latest_weekly_dataset.  After run_stage the
       lockfile is released
       :param theargs: should contain theargs.celppdir & other params
                       set via commandline
    """
    updatedtheargs = set_andor_create_latest_weekly_parameter(theargs)

    if updatedtheargs.latest_weekly is None:
        logger.info("No weekly dataset found in path " +
                    updatedtheargs.celppdir)
        return 0
    for stage_name in updatedtheargs.stage.split(','):
        logger.info("Starting " + stage_name + " stage")
        try:
            lock = _get_lock(updatedtheargs, stage_name)

            task_list = get_task_list_for_stage(updatedtheargs, stage_name)

            # run the stage
            exit_code = run_tasks(task_list)
            if exit_code is not 0:
                logger.error('Non zero exit code from task ' + stage_name +
                             'exiting')
                return exit_code
        finally:
            # release lock
            logger.debug('Releasing lock')
            lock.release()

    return 0


def run_tasks(task_list):
    """Runs a specific stage

       Runs the tasks in task_list
       :param task_list: list of tasks to run
    """
    if task_list is None:
        logger.error('Task list is None')
        return 3

    if len(task_list) == 0:
        logger.error('Task list is empty')
        return 2

    returnval = 0

    for task in task_list:
        logger.info("Running task " + task.get_name())
        try:
            task.run()
        except Exception as e:
            logger.exception("Error caught exception")
            if task.get_error() is None:
                task.set_error('Caught Exception running task: ' + e.message)

        logger.debug("Task " + task.get_name() + " has finished running " +
                     " with status " + task.get_status())
        if task.get_error() is not None:
            logger.error('Error running task ' + task.get_name() +
                         ' ' + task.get_error())
            returnval = 1

    return returnval


def get_task_list_for_stage(theargs, stage_name):
    """Factory method that generates a list of tasks for given stage

       Using stage_name get the list of tasks that need to
       be run.
       :param theargs: parameters set via commandline along with
                       ``theargs.latest_weekly`` which should be set to
                       to base directory where stages will be run
       :param stage_name:  Name of stage to run
    """
    if stage_name is None:
        raise NotImplementedError('stage_name is None')

    task_list = []

    logger.debug('Getting task list for ' + stage_name)

    if stage_name == CREATE_CHALLENGE:
        task_list.append(MakeBlastDBTask(theargs.latest_weekly, theargs))
        task_list.append(DataImportTask(theargs.latest_weekly, theargs))
        task_list.append(BlastNFilterTask(theargs.latest_weekly, theargs))
        task_list.append(ChallengeDataTask(theargs.latest_weekly, theargs))

    if stage_name == 'makedb':
        task_list.append(MakeBlastDBTask(theargs.latest_weekly, theargs))

    if stage_name == 'import':
        task_list.append(DataImportTask(theargs.latest_weekly, theargs))

    if stage_name == 'blast':
        task_list.append(BlastNFilterTask(theargs.latest_weekly, theargs))

    if stage_name == 'challengedata':
        task_list.append(ChallengeDataTask(theargs.latest_weekly, theargs))

    if stage_name == 'proteinligprep':
        task_list.append(ProteinLigPrepTask(theargs.latest_weekly, theargs))

    if stage_name == 'glide':
        task_list.append(GlideTask(theargs.latest_weekly, theargs))

    if stage_name == 'vina':
        task_list.append(AutoDockVinaTask(theargs.latest_weekly, theargs))

    if stage_name == CHIMERA_PREP:
        task_list.append(ChimeraProteinLigPrepTask(theargs.latest_weekly,
                                                   theargs))
    if stage_name == 'extsubmission':
        extfac = ExternalDataSubmissionFactory(theargs.latest_weekly, theargs)
        task_list.extend(extfac.get_external_data_submissions())

    if stage_name == 'evaluation':
        # use util function call to get all evaluation tasks
        # append them to the task_list
        eval_task_factory = EvaluationTaskFactory(theargs.latest_weekly,
                                                  theargs)
        task_list.extend(eval_task_factory.get_evaluation_tasks())

    if len(task_list) is 0:
        raise NotImplementedError(
            'uh oh no tasks for ' + stage_name + ' stage')

    return task_list


def _parse_arguments(desc, args):
    """Parses command line arguments using argparse.
    """
    parsed_arguments = D3RParameters()

    help_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=help_formatter)
    parser.add_argument("celppdir", help='Base celpp directory')
    parser.add_argument("--email", dest="email",
                        help='Comma delimited list of email addresses')
    parser.add_argument("--summaryemail", dest="summaryemail",
                        help='Comma delimited list of email addresses to'
                             'receive an email summarizing docking '
                             'evaluations DOES NOT DO ANYTHING YET')
    parser.add_argument("--createweekdir",
                        help='Create new celpp week directory before ' +
                             'running stages',
                        action="store_true")
    parser.add_argument("--customweekdir",
                        action='store_true',
                        help="Use directory set in celppdir instead of " +
                             "looking for latest weekdir.  NOTE: " +
                             "--createweekdir " +
                             "will create a dataset.week.# dir under celppdir")
    parser.add_argument("--stage", required=True,
                        help='Comma delimited list' +
                        ' of stages to run.  Valid STAGES = ' +
                        '{makedb, import, blast, challengedata,'
                        + CHIMERA_PREP +
                        ', proteinligprep, extsubmission, glide, '
                        'vina, evaluation, ' + CREATE_CHALLENGE + '} ')
    parser.add_argument("--blastnfilter", default='blastnfilter.py',
                        help='Path to BlastnFilter script '
                             '(default blastnfilter.py)')
    parser.add_argument("--blastnfiltertimeout", default=36000, type=int,
                        help='DOES NOT WORK WILL BE IGNORED. '
                             'Time in seconds script is allowed to run before'
                             'being killed (default 36000)')
    parser.add_argument("--postanalysis", default='postanalysis.py',
                        help='Path to PostAnalysis script '
                             '(default postanalysis.py)')
    parser.add_argument("--proteinligprep", default='proteinligprep.py',
                        help='Path to proteinligprep script (default '
                             'proteinligprep.py)')
    parser.add_argument("--genchallenge", default='genchallengedata.py',
                        help='Path to genchallengedata script '
                             '(default genchallengedata.py)')
    parser.add_argument("--chimeraprep", default='chimera_proteinligprep.py',
                        help='Path to chimera_proteinligprep script '
                             '(default chimera_proteinligprep.py)')
    parser.add_argument("--glide", default='glidedocking.py',
                        help='Path to glide docking script '
                             '(default glidedocking.py)')
    parser.add_argument("--vina", default='vinadocking.py',
                        help='Path to auto dock vina docking script '
                             '(default vinadocking.py)')
    parser.add_argument("--evaluation", default='evaluate.py',
                        help='Path to evaluation script (default evaluate.py)')
    parser.add_argument("--pdbdb", default='/data/pdb',
                        help='Path to PDB database files (default /data/pdb)')
    parser.add_argument("--compinchi",
                        default='http://ligand-expo.rcsb.org/' +
                        'dictionaries',
                        help='URL to download Components-inchi.ich' +
                             ' file for' +
                             'task stage dataimport '
                             '(default http://ligand-expo.rcsb.org/)')
    parser.add_argument("--pdbfileurl",
                        default='http://www.wwpdb.org/files',
                        help='URL to download ' +
                             DataImportTask.NONPOLYMER_TSV + ', ' +
                             DataImportTask.SEQUENCE_TSV + ', and ' +
                             DataImportTask.CRYSTALPH_TSV +
                             ' files for task stage dataimport '
                             '(default http://www.wwpdb.org/files)')
    parser.add_argument("--makeblastdb", default='makeblastdb',
                        help='Path to NCBI Blast makeblastdb program '
                             'ie /usr/bin/makeblastdb (default makeblastdb)')
    parser.add_argument("--pdbsequrl",
                        default='ftp://ftp.rcsb.org/pub/pdb/derived_data/'
                                'pdb_seqres.txt.gz',
                        help='ftp url to download rcsb sequences file '
                             '(default '
                             'ftp://ftp.rcsb.org/pub/pdb/derived_data/'
                             'pdb_seqres.txt.gz')
    parser.add_argument("--rdkitpython",
                        default='',
                        help='Path to directory containing python with '
                             'new rdkit installed')
    parser.add_argument("--log", dest="loglevel", choices=['DEBUG',
                        'INFO', DEFAULT_LOG_LEVEL, 'ERROR', 'CRITICAL'],
                        help="Set the logging level (default " +
                             DEFAULT_LOG_LEVEL + ")",
                        default=DEFAULT_LOG_LEVEL)
    parser.add_argument('--smtp', dest='smtp', help='Sets smtpserver to use '
                                                    '(default localhost)',
                        default='localhost')
    parser.add_argument('--smtpport', dest='smtpport',
                        help='Sets smtp server port (default 25)',
                        default='25')
    parser.add_argument('--skipimportwait', default=False, action='store_true',
                        help='Normally the import stage will wait if any of '
                             'the tsv files have not been updated since the'
                             'start of the current celpp week.  Setting this '
                             'flag bypasses this delay and downloads the'
                             'files')
    parser.add_argument('--importsleep', default=600, type=int,
                        help='Number of seconds to wait before re-checking '
                             'tsv files to see if they have been updated in '
                             'the import stage (default 600)')
    parser.add_argument('--importretry', default=60, type=int,
                        help='Number of times import stage should check if '
                             'tsv files have been updated (default 60, or 10 '
                             'hours with default --importsleep set to 600 '
                             'seconds)')
    parser.add_argument('--ftpconfig', dest='ftpconfig', help='File containing'
                        ' configuration to connect to ftp server.  If set,'
                        ' data from stages run during this invocation will be'
                        ' uploaded after the stage completes.  Format is same'
                        ' as ncftp config files with two added fields (' +
                        FtpFileTransfer.PATH + ', ' +
                        FtpFileTransfer.CHALLENGEPATH + ', ' +
                        FtpFileTransfer.SUBMISSIONPATH + ') ' +
                        'SEE: description of challengedata '
                        'stage above for example file')
    parser.add_argument('--replyto', dest='replytoaddress', default=None,
                        help='Email address to set as reply-to for emails to '
                             'external submitters')
    parser.add_argument('--version', action='version',
                        version=('%(prog)s ' + d3r.__version__))
    return parser.parse_args(args, namespace=parsed_arguments)


def main(args):
    """Main entry into celpprunner
    :param args: should be set to sys.argv which is a list of arguments
                 starting with script name as the first argument
    """
    blasttask = BlastNFilterTask('', p)
    dataimport = DataImportTask('', p)
    challenge = ChallengeDataTask('', p)
    glide = GlideTask('', p)
    makedb = MakeBlastDBTask('', p)
    prot = ProteinLigPrepTask('', p)
    vina = AutoDockVinaTask('', p)
    chimeraprep = ChimeraProteinLigPrepTask('', p)
    desc = """
              Version {version}

              Runs the 9 stages (makedb, import, blast, challengedata,
              proteinligprep, {chimeraprep}, extsubmission, glide, vina, &
              evaluation) of CELPP processing pipeline
              (http://www.drugdesigndata.org)

              CELPP processing pipeline relies on a set of directories
              with specific structure. The pipeline runs a set of stages
              Each stage has a numerical value and a name. The numerical
              value denotes order and the stage name identifies separate
              tasks to run in the stage.

              The filesystem structure of the stage is:

              stage.<stage number>.<task name>

              The stage(s) run are defined via the required --stage flag.

              To run multiple stages serially just pass a comma delimited
              list to the --stage flag. Example: --stage import,blast

              NOTE:  When running multiple stages serially the program will
                     not run subsequent stages if a task in a stage fails.
                     Also note order matters, ie putting blast,import will
                     cause celpprunner.py to run blast stage first.

              This program drops a pid lockfile
              (celpprunner.<stage>.lockpid) in celppdir to prevent duplicate
              invocation.

              When run, this program will examine the stage and see
              if work can be done.  If stage is complete or previous
              steps have not completed, the program will exit silently.
              If previous steps have failed or current stage already
              exists in an error or uncomplete state then program will
              report the error via email using addresses set in --email
              flag. Errors will also be reported via stderr/stdout.
              The program will also exit with nonzero exit code.

              This program utilizes simple token files to denote stage
              completion.  If within the stage directory there is a:

              '{complete}' file - then stage is done and no other
                                checking is done.

              'error' file - then stage failed.

              'start' file - then stage is running.

              Notification of stage start and end will be sent to
              addresses set via --email flag.

              Unless --customweekdir is set, this program will
              examine the 'celppdir' (last argument passed on
              commandline) to find the latest directory with this path:

              <year>/dataset.week.#

              The program will find the latest <year> and within
              that year the dataset.week.# with highest #.  The output
              directories created will be put within this directory.

              Setting --customweekdir will cause program to use 'celppdir'
              path.

              Setting the --createweekdir flag will instruct this
              program to create a new directory for the current
              celpp week/year before running any stage processing.

              NOTE: CELPP weeks start on Friday and end on Thursday
                    and week # follows ISO8601 rules so week numbers
                    at the end and start of the year are a bit
                    wonky.

              Breakdown of behavior of program is defined by
              value passed with --stage flag:

              If --stage '{createchallenge}'

              This is NOT a stage, but has the same effect as
              calling --stage makedb,import,blast,challengedata
              The four stages that need to run to generate the challenge
              data package.

              If --stage 'makedb'

              In this stage the file {pdb_seqres} is downloaded from
              an ftp site set by --pdbsequrl.
              This file is then gunzipped and NCBI makeblastdb
              (set by --makeblastdb) is run on it to create a blast
              database.  The files are stored in {makeblastdb_dirname}

              If --stage 'import'

              In this stage 4 files are downloaded from urls specified
              by --compinchi and --pdbfileurl flags on the commandline
              into {dataimport_dirname} directory.

              The tsv files are (--pdbfileurl flag sets url to
              download these files from):

              {nonpolymer_tsv}
              {sequence_tsv}
              {crystal_tsv}

              The Components ich file is (--compinchi flag sets base url to
              download this file from):

              {compinchi_ich}

              This stage will just wait and retry if any of the tsv files
              have NOT been updated since the start of the current
              celpp week as determined by a HEAD request. To bypass
              this delay add --skipimportwait flag.  --importsleep denotes
              the time to wait before re-examining the update time of the
              tsv files and --importretry sets number of times to retry
              before giving up.

              If --stage 'blast'

              Verifies {dataimport_dirname} exists and has '{complete}'
              file.  Also verifies {makeblastdb_dirname} exists and has
              '{complete}' file.  If both conditions are met then the
              'blast' stage is run which invokes script set by
              --blastnfilter flag and output stored in
              {blast_dirname}.
              Requires --pdbdb to be set to a directory with valid PDB
              database files.

              Note: --blastnfilter script is killed after time set with
              --blastnfiltertimeout flag.


              If --stage 'challengedata'

              Verifies {blast_dirname} exists and has '{complete}'
              file.  If complete, this stage runs which invokes program
              set in --genchallenge flag to create a challenge dataset
              file.  The --pdbdb flag must also be set when calling this
              stage. If --ftpconfig is set with {challengepath} field then
              this stage will also upload the challenge dataset tarfile
              to the ftp server with path set by {challengepath}.  The
              code will also upload a {latest_txt} file containing name
              of the tarfile to the same destination overwriting any
              {latest_txt} file that already exists.

              Example file for --ftpconfig:

              {host} some.ftp.com
              {user} bob
              {passn} mypass
              {path} /celpp
              {challengepath} /challenge
              {submissionpath} /submissions


              If --stage '{chimeraprep}'

              Verifies {challenge_dirname} exists and has '{complete}'
              file.  If complete, this stage runs which invokes program
              set in --chimeraprep flag to prepare pdb and inchi files
              storing output in {chimeraprep_dirname}.  --pdbdb flag
              must also be set when calling this stage.

              If --stage 'proteinligprep'

              Verifies {challenge_dirname} exists and has '{complete}'
              file.  If complete, this stage runs which invokes program
              set in --proteinligprep flag to prepare pdb and inchi files
              storing output in {proteinligprep_dirname}.  --pdbdb flag
              must also be set when calling this stage.

              If --stage 'extsubmission'

              Connects to server specified by --ftpconfig and downloads
              external docking submissions from {submissionpath} on remote
              server.

              Submissions should be named:

              celpp_weekXX_YYYY_dockedresults_ZZZZ.tar.gz as documented here:

              https://github.com/drugdata/d3r/wiki/Proposed-challenge-docked\
              -results-file-structure

              For each submission a directory named stage.X.ZZZZ.extsubmission
              will be created and uncompressed contents of package will be
              stored in that directory.  If data does not conform properly
              'error' file will be placed in directory denoting failure

              If --stage 'glide'

              Verifies {proteinligprep_dirname} exists and has a '{complete}'
              file within it.  If complete, this stage runs which invokes
              program set in --glide flag to perform docking via glide
              storing output in {glide_dirname}

              If --stage 'vina'

              Verifies {proteinligprep_dirname} exists and has a '{complete}'
              file within it.  If complete, this stage runs which invokes
              program set in --vina flag to perform docking via AutoDock Vina
              storing output in {vina_dirname}

              If --stage 'evaluation'

              Finds all stage.{dockstage}.<algo> directories with '{complete}'
              files in them which do not end in name '{webdata}' and runs
              script set via --evaluation parameter storing the result of
              the script into stage.{evalstage}.<algo>.evaluation. --pdbdb flag
              must also be set when calling this stage.


              """.format(makeblastdb_dirname=makedb.get_dir_name(),
                         dataimport_dirname=dataimport.get_dir_name(),
                         blast_dirname=blasttask.get_dir_name(),
                         challenge_dirname=challenge.get_dir_name(),
                         createchallenge=CREATE_CHALLENGE,
                         proteinligprep_dirname=prot.get_dir_name(),
                         glide_dirname=glide.get_dir_name(),
                         vina_dirname=vina.get_dir_name(),
                         dockstage=str(glide.get_stage()),
                         evalstage=str(glide.get_stage() + 1),
                         complete=blasttask.COMPLETE_FILE,
                         chimeraprep_dirname=chimeraprep.get_dir_name(),
                         chimeraprep=CHIMERA_PREP,
                         compinchi_ich=DataImportTask.COMPINCHI_ICH,
                         pdb_seqres=MakeBlastDBTask.PDB_SEQRES_TXT_GZ,
                         nonpolymer_tsv=DataImportTask.NONPOLYMER_TSV,
                         sequence_tsv=DataImportTask.SEQUENCE_TSV,
                         crystal_tsv=DataImportTask.CRYSTALPH_TSV,
                         webdata=EvaluationTaskFactory.WEB_DATA_SUFFIX,
                         latest_txt=ChallengeDataTask.LATEST_TXT,
                         host=FtpFileTransfer.HOST,
                         user=FtpFileTransfer.USER,
                         passn=FtpFileTransfer.PASS,
                         path=FtpFileTransfer.PATH,
                         challengepath=FtpFileTransfer.CHALLENGEPATH,
                         submissionpath=FtpFileTransfer.SUBMISSIONPATH,
                         version=d3r.__version__)

    theargs = _parse_arguments(desc, args[1:])
    theargs.program = args[0]
    theargs.version = d3r.__version__

    util.setup_logging(theargs)

    try:
        return run_stages(theargs)
    except Exception:
        logger.exception("Error caught exception")
        return 2


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
