#! /usr/bin/env python

import sys
import os
import argparse
import psutil
import logging

from d3r import util
from d3r.task import D3RParameters
from d3r.task import BlastNFilterTask
from lockfile.pidlockfile import PIDLockFile

# create logger
logger = logging.getLogger('d3r.celpprunner')
LOG_FORMAT = "%(asctime)-15s %(levelname)s %(name)s %(message)s"


def _get_lock(theargs):
    """Create lock file to prevent this process from running on same data

       This uses ``PIDLockFile`` to create a pid lock file in celppdir
       directory named celprunner.<stage>.lockpid
       If pid exists it is assumed the lock is held otherwise lock
       is broken and recreated

       :param theargs: return value from argparse and should contain
                       theargs.stage which denotes stage of processing
                       and theargs.celppdir should be set to path
       :return: ``PIDLockFile`` upon success
       :raises: LockException
       """
    mylockfile = os.path.join(theargs.celppdir, "celpprunner." +
                              theargs.stage + ".lockpid")
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


def _setup_logging(theargs):
    """Sets up the logging for application
       """
    theargs.logformat = LOG_FORMAT
    logger.setLevel(theargs.loglevel)
    logging.basicConfig(format=theargs.logformat)
    logging.getLogger('d3r.task').setLevel(theargs.loglevel)


def run_stage(theargs):
    theargs.latest_weekly = util.find_latest_weekly_dataset(theargs.celppdir)

    if theargs.latest_weekly is None:
        logger.info("No weekly dataset found in path " +
                    theargs.celppdir)
        return 0

    logger.info("Starting " + theargs.stage + " stage")

    # perform processing
    if theargs.stage == 'blast':
        task = BlastNFilterTask(theargs.latest_weekly, theargs)

    if theargs.stage == 'dock':
        raise NotImplementedError('uh oh dock is not implemented yet')

    if theargs.stage == 'score':
        raise NotImplementedError('uh oh score is not implemented yet')

    logger.info("Running task " + task.get_name())
    task.run()
    logger.debug("Task " + task.get_name() + " has finished running " +
                 " with status " + task.get_status())
    if task.get_error() is not None:
        logger.error('Error running task ' + task.get_name() +
                     ' ' + task.get_error())
        return 1
    return 0


def _parse_arguments(desc, args):
    """Parses command line arguments
       """
    parsed_arguments = D3RParameters()

    help_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=help_formatter)
    parser.add_argument("celppdir", help='Base celpp directory')
    parser.add_argument("--blastdir", help='Parent directory of ' +
                        ' blastdb.  There should exist a "current" ' +
                        ' symlink or directory that contains the db.')
    parser.add_argument("--email", dest="email",
                        help='Comma delimited list of email addresses')

    parser.add_argument("--stage", choices=['blast', 'dock', 'score'],
                        required=True, help='Stage to run blast = ' +
                        'blastnfilter (2), dock = fred & other ' +
                        'docking algorithms (3), ' +
                        'score = scoring (4)')
    parser.add_argument("--blastnfilter", required=True,
                        help='Path to BlastnFilter script')
    parser.add_argument("--log", dest="loglevel", choices=['DEBUG',
                        'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help="Set the logging level",
                        default='WARNING')
    parser.add_argument('--smtp', dest='smtp', help='Sets smtpserver to use',
                        default='localhost')
    parser.add_argument('--smtpport', dest='smtpport',
                        help='Sets smtp server port', default='25')

    return parser.parse_args(args, namespace=parsed_arguments)


def main():

    desc = """
              Runs last 3 stages (blast, dock, & score) of CELPP
              processing pipeline (http://www.drugdesigndata.org)

              CELPP processing pipeline relies on a set of directories
              with specific structure. The pipeline runs a set of stages
              Each stage has a numerical value and a name. The numerical
              value denotes order and the stage name identifies separate
              tasks to run in the stage.

              The filesystem structure of the stage is:

              stage.<stage number>.<task name>

              Only 1 stage is run per invocation of this program and the
              stage to be run is defined via the required --stage flag.

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

              'complete' file - then stage is done and no other
                                checking is done.

              'error' file - then stage failed.

              'start' file - then stage is running.

              Notification of stage start and end will be sent to
              addresses set via --email flag.

              Regardless of the stage specified, this program will
              examine the 'celppdir' (last argument passed on
              commandline) to find the latest directory with this path:
              <year>/dataset.week.#
              The program will find the latest <year> and within
              that year the dataset.week.# with highest #.  The output
              directories created will be put within this directory.

              Breakdown of behavior of program is defined by
              value passed with --stage flag:

              If --stage 'blast'

              Verifies stage.1.dataimport exists and has 'complete'
              file.  Also the --blastdir path must exist and within a
              'current' symlink/directory must exist and within that a
              'complete' file must also reside. If both conditions
              are met then the 'blast' stage is run and output stored
              in stage.2.blastnfilter

              If --stage 'dock'

              Verifies stage2.blastnfilter exists and has a 'complete'
              file within it.  If complete, this program will run fred
              docking and store output in stage.3.fred.  As new
              algorithms are incorporated additional stage.3.<algo> will
              be created and run.

              If --stage 'score'

              Finds all stage.3.<algo> directories with 'complete' files
              in them and invokes appropriate scoring algorithm storing
              results in stage.4.<algo>.scoring.
              """

    theargs = _parse_arguments(desc, sys.argv[1:])
    theargs.program = sys.argv[0]
    try:
        if os.path.basename(theargs.blastdir) is 'current':
            theargs.blastdir = os.path.dirname(theargs.blastdir)
    except AttributeError:
        pass

    _setup_logging(theargs)

    try:
        # get the lock
        lock = _get_lock(theargs)

        # run the stage
        sys.exit(run_stage(theargs))

    except Exception:
        logger.exception("Error caught exception")
        sys.exit(2)
    finally:
        # release lock
        logger.debug('Releasing lock')
        lock.release()


if __name__ == '__main__':
    main()
