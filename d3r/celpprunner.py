#! /usr/bin/env python

import os
import argparse
import psutil
import logging

from d3r.task import D3RTask
from lockfile.pidlockfile import PIDLockFile

# create logger
logger = logging.getLogger('d3r.celpprunner')
LOG_FORMAT = "%(asctime)-15s %(levelname)s %(name)s %(message)s"

def get_lock(theargs):
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
   mylockfile = os.path.join(theargs.celppdir,"celpprunner."+
                                              theargs.stage+".lockpid")
   logger.debug("Looking for lock file: "+mylockfile)
   lock = PIDLockFile(mylockfile,timeout=10)
   
   if lock.i_am_locking():
       logger.debug("My process id"+str(lock.read_pid())+
                    " had the lock so I am breaking")
       lock.break_lock()
       lock.acquire(timeout=10)
       return lock

   if lock.is_locked():
       logger.debug("Lock file exists checking pid")
       if psutil.pid_exists(lock.read_pid()):
           raise Exception("celpprunner with pid "+str(lock.read_pid())+
                           " is running")

   lock.break_lock()
   logger.info("Acquiring lock")
   lock.acquire(timeout=10)
   return lock

def setup_logging(theargs):
   """Sets up the logging for application
      """
   theargs.logFormat = LOG_FORMAT
   logger.setLevel(theargs.logLevel)

   logging.basicConfig(format=theargs.logFormat)


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

   parser = argparse.ArgumentParser(description=desc)
   parser.add_argument("celppdir",help='Base celpp directory')
   parser.add_argument("--stage", choices=['blast','dock','score'],
                       required=True,help='Stage to run blast = '+
     'blastnfilter (2), dock = fred & other docking algorithms (3), '+
     'score = scoring (4)')  
   parser.add_argument("--log", dest="logLevel", choices=['DEBUG', 
                       'INFO', 'WARNING', 'ERROR', 'CRITICAL'], 
                       help="Set the logging level",
                       default='WARNING')
 
   theargs = parser.parse_args()

   setup_logging(theargs)
 
   # get the lock
   lock = get_lock(theargs)

     
   # release lock
   lock.release()

#  print "Hello World my process id is: ",mypid
#  task = D3RTask()
#  print task.get_name()

if __name__ == '__main__':
  main()

