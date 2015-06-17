#! /usr/bin/env python

import sys
import os
import argparse
import psutil
import logging

import d3r.task
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
    theargs.logFormat = LOG_FORMAT
    logger.setLevel(theargs.logLevel)
    logging.basicConfig(format=theargs.logFormat)


def _parse_arguments(desc, args):
    """Parses command line arguments
       """
    parsed_arguments = D3RParameters()

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("celppdir", help='Base celpp directory')
    parser.add_argument("--blastdir", help='Directory containing blastdb')
    parser.add_argument("--email",
                        help='Comma delimited list of email addresses')

    parser.add_argument("--stage", choices=['blast', 'dock', 'score'],
                        required=True, help='Stage to run blast = ' +
                        'blastnfilter (2), dock = fred & other ' +
                        'docking algorithms (3), ' +
                        'score = scoring (4)')
    parser.add_argument("--log", dest="logLevel", choices=['DEBUG',
                        'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help="Set the logging level",
                        default='WARNING')

    return parser.parse_args(args, namespace=parsed_arguments)


def main():

    desc = """Runs last 3 stages of CELPP processing pipeline (blast,
              docking, and scoring).  This tool will examine the
              celppdir to find the latest weekly download of data from
              wwPDB which should be under\n
              <year>/dataset.week.#/stage.1.dataimport path.\n
              The tool then verifies the stage specified by can be run
              and performs the operation to prevent duplicate invocation
              a token file named celprunner.# is dropped in the celppdir
              which contains the pid of the process.  This is checked
              upon startup to prevent duplicate invocation.
              """

    theargs = _parse_arguments(desc, sys.argv[1:])

    _setup_logging(theargs)

    # get the lock
    lock = _get_lock(theargs)

    latestWeekly = d3r.task.find_latest_weekly_dataset(theargs.celppdir)

    if latestWeekly is None:
        logger.debug("No weekly dataset found")
        return

    # perform processing
    if theargs.stage == 'blast':
        print "Blast stage"
        task = BlastNFilterTask(theargs, latestWeekly)

    if theargs.stage == 'dock':
        print "dock stage"
        return

    if theargs.stage == 'score':
        print "score stage"
        return

    if not task.can_run():
        logger.debug("Task " + task.get_name() + " cannot run ")

    logger.debug("Running task " + task.get_name())
    task.run()
    logger.debug("Task " + task.get_name() + " has finished running " +
                 " with status " + task.get_status())

    # release lock
    lock.release()


if __name__ == '__main__':
    main()
