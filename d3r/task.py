# -*- coding: utf-8 -*-

import os
import logging
import re

logger = logging.getLogger(__name__)


class D3RParameters(object):
    pass


class D3RTask:
    """Represents a base Task that can be run.

    This is a base class from which other classes that actual do work
    can be derived from.  It provides internal attributes for name,
    stage, and status.  As well as a run() function.

    """

    START_STATUS = "start"
    COMPLETE_STATUS = "complete"
    UNKNOWN_STATUS = "unknown"
    NOTFOUND_STATUS = "notfound"
    ERROR_STATUS = "error"

    def __init__(self, path, args):
        self._path = path
        self._name = None        
        self._stage = None
        self._status = D3RTask.UNKNOWN_STATUS
        self._error = None
        self._args = args

    def get_path(self):
        return self._path

    def set_path(self, path):
        self._path = path

    def get_name(self):
        return self._name

    def set_name(self, name):
        self._name = name

    def get_stage(self):
        return self._stage

    def set_stage(self, stage):
        self._stage = stage

    def get_status(self):
        return self._status

    def set_status(self, status):
        self._status = status

    def get_error(self):
        return self._error

    def set_error(self, error):
        self._error = error

    def start(self):
        print "Task ", self._name, " is starting"

    def end(self):
        print "Task ", self._name, " has finished"

    def run(self):
        print "Running ", self._name

    def get_dir_name(self):
        """Gets directory name for task

           """
        return 'stage' + '.' + self._stage + '.' + self._name

    def update_status_from_filesystem(self):
        """Updates status by querying filesystem.

           Sets status based on contents of path on filesystem.
           If path does not exist then status is NOTFOUND_STATUS.
           If complete file exists under path then status is COMPLETE_STATUS
           If error file exists under path then status is ERROR_STATUS
           If start file exists under path then status is START_STATUS
           else status is UNKNOWN_STATUS       
           
           """
        pass


    def create_dir(self):
        """Creates directory for task.
            
           Directory will be named stage.<stage>.<name>
           and located under get_path()
           """
        return os.mkdir(os.path.join(self.get_path(), self.get_dir_name()))

class DataImportTask(D3RTask):
    """Represents DataImport Task

       """

    def __init__(self, path, args):
        D3RTask.__init__(self, path, args)
        self.set_name('dataimport')
        self.set_stage(1)
        self.set_status(D3RTask.UNKNOWN_STATUS)


class MakeBlastDBTask(D3RTask):
    """Represents creation of BlastDB

       """

    def __init__(self, path, args):
        D3RTask.__init__(self, path, args)
        self.set_name('makeblastdb')
        self.set_stage(0)
        self.set_status(D3RTask.UNKNOWN_STATUS)

    def update_status_from_filesystem(self):
        """Updates status by querying filesystem.

           Sets status based on contents of path on filesystem.
           If path does not exist then status is NOTFOUND_STATUS.
           If complete file exists under path then status is COMPLETE_STATUS
           If error file exists under path then status is ERROR_STATUS
           If start file exists under path then status is START_STATUS
           else status is UNKNOWN_STATUS       
           
           """
        pass
     

class BlastNFilterTask(D3RTask):
    """Performs Blast and filter of sequences

    """

    def __init__(self, path, args):
        D3RTask.__init__(self, path, args)
        self.set_name('blastnfilter')
        self.set_stage(2)
        self.set_status(D3RTask.UNKNOWN_STATUS)

    def can_run(self):
        
        # check blast 
        makeblastdb = MakeBlastDBTask(self._path,self._args)
        makeblastdb.update_status_from_filesystem()
        if self.get_status() != D3RTask.COMPLETE_STATUS:
            return False

        # check data import
        dataImport = DataImportTask(self._path,self._args)
        dataImport.update_status_from_filesystem()
        if self.get_status() != D3RTask.COMPLETE_STATUS:
            return False

        # check blast is not complete and does not exist

        self.update_status_from_filesystem()
        if self.get_status() == D3RTask.COMPLETE_STATUS:
            logger.debug("No work needed " + self.get_name() +
                         " task is complete")
            return False

        if self.get_status() != D3R.NOTFOUND_STATUS:
            logger.debug("Task was already attempted, but there was a problem")
            return False 

        return True

    def run(self):
        """Runs blastnfilter task after verifying dataimport was good

           First checks that stage.1.dataimport has correct data and 
           that blast db is valid and that blastnfilter was not already
           run in any form.  If not method just returns doing nothing.
           Otherwise method invokes D3RTask.start then this method
           creates a directory and invokes blastnfilter script.  Upon
           completion results are analyzed and success or error status
           is set appropriately and D3RTask.end is invoked
           """
        logger.debug("In BlastNFilter Task")

        if not self.can_run():
            return

        self.start()
        
        if not self.create_dir():
            logger.error("Unable to create directory")
            self.set_error("Unable to create directory")
            self.end() 
            return

        # Run the blastnfilter
        print "blastnfilter <> <>"

        # assess the result
        self.end() 


def find_latest_year(celppdir):
    """Given a directory find the latest year

       The latest year will be a folder of 4 digits
       and have the highest value ie 2015
       :return: Directory or None if none is found
       :raises: Exception: If celppdir is not a directory
       """

    if not os.path.isdir(celppdir):
        raise Exception(celppdir+" is not a directory")

    dirPattern = re.compile("^\d\d\d\d$")

    latestYear = -1
    latestEntry = None
    fullPath = None
    for entry in os.listdir(celppdir):
        if re.match(dirPattern, entry):
            entryYear = int(entry)

            if entryYear > latestYear:
                fullPath = os.path.join(celppdir, entry)
                if os.path.isdir(fullPath):
                    latestYear = entryYear
                    latestEntry = fullPath

    if latestYear == -1:
        return

    return latestEntry


def find_latest_weekly_dataset(celppdir):
    """Given a directory find the latest dataset

       This method looks in directory passed in for
       paths with dataset.week.# format and returns
       the path with highest # on it.  If none are
       found then None is returned.
       :return: Directory upon success or None if no directory found
       :raises: Exception: If celppdir is not a directory
       """

    latestYear = find_latest_year(celppdir)

    if latestYear is None:
        return

    celppdirWithYear = os.path.join(celppdir, latestYear)

    dirPattern = re.compile("^dataset.week.\d+$")

    latestEntry = None
    latestWeekNo = -1
    fullPath = None
    for entry in os.listdir(celppdirWithYear):
        if re.match(dirPattern, entry):
            weekNo = re.sub("dataset.week.", "", entry)
            if weekNo > latestWeekNo:
                fullPath = os.path.join(celppdirWithYear, entry)
                if os.path.isdir(fullPath):
                    latestWeekNo = weekNo
                    latestEntry = fullPath

    return latestEntry
