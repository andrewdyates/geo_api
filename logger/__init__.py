#!/usr/bin/python
import stderr_logger

# set context sensitive logger
if "SERVER" in globals() and SERVER:
    Log = server_logger.Log
else:
    Log = stderr_logger.Log
