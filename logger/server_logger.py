#!/usr/bin/python
"""Logger which writes to log file.

For help, see:

"Logging Cookbook"
Python v2.7.2 documentation: Python HOWTOs.
[http://docs.python.org/howto/logging-cookbook.html]
"""
import logging


Log = logging.getLogger("server_logger")
Log.setLevel(logging.DEBUG)

# StreamHandler() uses sys.stderr by default.
# http://docs.python.org/library/logging.handlers.html
LOGFILE = 'log' # XXX: THIS NEEDS TO BE SET DYNAMICALLY
fh = logging.FileHandler(LOGFILE)
fh.setLevel(logging.DEBUG)
f = logging.Formatter("%(levelname)s %(asctime)s %(module)s.%(funcName)s %(lineno)d %(message)s")
fh.setFormatter(f)
Log.addHandler(fh)

def console(s):
  pass

Log.console = console
