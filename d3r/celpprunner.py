#! /usr/bin/env python

import os

from d3r.task import D3RTask

def main():
  mypid = os.getpid()

  print "Hello World my process id is: ",mypid

  task = D3RTask()
  print task.get_name()

if __name__ == '__main__':
  main()

