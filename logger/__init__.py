#!/usr/bin/python
import os
import stderr_logger

# Set context sensitive logger.
if 'ENV' in os.environ and os.environ['ENV'] == "SERVER":
    Log = server_logger.Log
else:
    Log = stderr_logger.Log
