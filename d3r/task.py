# -*- coding: utf-8 -*-

import os
import logging
import re

from os import listdir

logger = logging.getLogger(__name__)

class D3RParameters(object):
    pass

class D3RTask:
   """Represents a base Task that can be run.
   
   This is a base class from which other classes that actual do work
   can be derived from.  It provides internal attributes for name,
   stage, and status.  As well as a run() function.  

   """

   def __init__(self,args):
       self._name = "hi"
       self._stage = 1
       self._status = "unknown"
       self._error = ""
       self._args = args

   def get_name(self):
       return self._name

   def set_name(self,name):
       self._name = name

   def get_stage(self):
       return self._stage

   def set_stage(self,stage):
       self._stage = stage

   def get_status(self):
       return self._status

   def set_status(self,status):
       self._status = status

   def get_error(self):
       return self._error

   def set_error(self,error):
       self._error = error

   def start(self):
       print "Task ",self._name," is starting"

   def end(self):
       print "Task ",self._name," has finished"

   def run(self):
       print "Running ",self._name

class BlastNFilterTask(D3RTask):
   """Performs Blast and filter of sequences


   """

   def __init__(self,args):
       D3RTask.__init__(self,args)

   def run(self):
       print "BlastNFilterTask"
       logger.debug("In BlastNFilter Task")


def _find_latest_year(celppdir):
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
        if re.match(dirPattern,entry):
            entryYear = int(entry)
            
            if entryYear > latestYear:
                fullPath = os.path.join(celppdir,entry)
                if os.path.isdir(fullPath):
                    latestYear = entryYear
                    latestEntry = fullPath
                
    if latestYear == -1:
        return

    return latestEntry
    

def _find_latest_weekly_dataset(celppdir):
    """Given a directory find the latest dataset
  
       This method looks in directory passed in for
       paths with dataset.week.# format and returns
       the path with highest # on it.  If none are
       found then None is returned.
       :return: Directory upon success or None if no directory found
       :raises: Exception: If celppdir is not a directory
       """

    latestYear = _find_latest_year(celppdir)

    if latestYear is None:
        return

    celppdirWithYear = os.path.join(celppdir,latestYear)

    dirPattern = re.compile("^dataset.week.\d+$")

    latestEntry = None
    latestWeekNo = -1
    for entry in os.listdir(celppdirWithYear):
        if re.match(dirPattern,entry):
           weekNo = re.sub("dataset.week.","",entry)
           if weekNo > latestWeekNo:
               latestWeekNo = weekNo
               latestEntry = entry

    if latestWeekNo == -1:
        return
  
    return os.path.join(celppdirWithYear,latestEntry)
