import sys
import os
import time

# HANDLE GLOBAL ENVIRONMENT AUTOMATICALLY AND BY DEFAULT
# --------------------
# Verify that environment variables are set. 
# If NONE of them are set, we may assume that the user would like to 
#   use default values: the local, current directory environment.
if ("ENV" not in os.environ) and ("CACHE_DIR" not in os.environ) and \
  ("TMP_DIR" not in os.environ):
  print "Warning: geo_api environment is not configured. Using local directory..."
  os.environ["ENV"] = "LOCAL"
  os.environ["CACHE_DIR"] = ""
  os.environ["TMP_DIR"] = ""

from __init__ import *


gse_id = "GSE25935"
gpl_id = None
g = GSE(gse_id, platform_id=gpl_id)

for name, sample in g.samples.items():
  row = {'id': name}
  row.update(dict([s.split(': ') for s in sample.attr['characteristics_ch1']]))
  for k,v in row.items():
    try:
      q = float(v)
    except ValueError:
      continue
    else:
      if "." in v:
        row[k] = q
      else:
        row[k] = int(v)
  print row


