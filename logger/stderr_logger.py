#!/usr/bin/python
"""A simple, noisy logging handle which prints to STDERR.

Intended for local development.

For help, see:

"Logging Cookbook"
Python v2.7.2 documentation: Python HOWTOs.
[http://docs.python.org/howto/logging-cookbook.html]
"""
import os
import logging

Log = logging.getLogger("stderr_logger")
Log.setLevel(logging.DEBUG)

# StreamHandler() uses sys.stderr by default.
# http://docs.python.org/library/logging.handlers.html
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
f = logging.Formatter("%(levelname)s %(asctime)s %(module)s.%(funcName)s %(lineno)d %(message)s")
ch.setFormatter(f)
Log.addHandler(ch)

# Create file handler which logs even debug messages.
# Use "debug" for big file dumps
home_dir = os.environ["HOME"]
fh = logging.FileHandler("%s/%s" % (home_dir, 'python_logger.log'))
fh.setLevel(logging.DEBUG)
Log.addHandler(fh)

import sys
def console(s):
  sys.stdout.write(s + '\r')
  sys.stdout.flush()
  
Log.console = console
