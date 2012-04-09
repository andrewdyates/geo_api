#!/usr/bin/python
import os

# Set context sensitive logger.
if 'ENV' in os.environ and os.environ['ENV'] == "SERVER":
    import server_logger
    Log = server_logger.Log
else:
    import stderr_logger
    Log = stderr_logger.Log
