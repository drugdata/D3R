# -*- coding: utf-8 -*-

class D3RTask:
   """Represents a base Task that can be run.
   
   This is a base class from which other classes that actual do work
   can be derived from.  It provides internal attributes for name,
   stage, and status.  As well as a run() function.  

   """

   def __init__(self):
       self._name = "hi"
       self._stage = 1
       self._status = "unknown"
       self._error = ""

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
   def __init__(self):
       D3RTask.__init__(self)

   def run(self):
       print "Hello me"
