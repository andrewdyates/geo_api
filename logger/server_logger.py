#!/usr/bin/python
"""Logger which writes to log file.

For help, see:

"Logging Cookbook"
Python v2.7.2 documentation: Python HOWTOs.
[http://docs.python.org/howto/logging-cookbook.html]
"""
import logging
import os

# Get log directory from environment variable.
if 'LOGFILE' not in os.environ:
  raise Exception, "Set os.environ variable 'LOGFILE' to full path to log file."
else:
  log_file = os.environ['LOGFILE']

Log = logging.getLogger(log_file)
Log.setLevel(logging.DEBUG)

# See: http://docs.python.org/library/logging.handlers.html
fh = logging.FileHandler(LOGFILE)
fh.setLevel(logging.DEBUG)
f = logging.Formatter("%(levelname)s %(asctime)s %(module)s.%(funcName)s %(lineno)d %(message)s")
fh.setFormatter(f)
Log.addHandler(fh)

def console(s):
  pass

Log.console = console
