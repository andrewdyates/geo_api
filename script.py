#!/usr/bin/python
"""Download a GSE Study in a clean format. Writes to local directory by default.
Packaged as an easy script wrapper that doesn't require the local environment
variables to be set and intelligently handles substudies.

Downloads each substudy as a separate file. Writes log about results.

SAMPLE USE:
$ python script.py GSE15745 2> api_log.txt
"""
USE_MSG = """USE: python script.py GSE_ID [options]

OPTIONS:
  GPL_ID=str: pseudo substudy GPL ID
  out_dir=str: path to an output directory
"""

import sys
import os
import time

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


def report(msg, fp):
  """Both write a message to a log file and to the console."""
  fp.write(msg + "\n")
  print msg


def main(gse_id, gpl_id=None, out_dir=""):

  # Verify that out_dir exists, and if not, create it.
  if not (os.path.exists(out_dir) and os.path.isdir(out_dir)):
    os.makedirs(out_dir)
  # Open log, write job description
  fp_log = open(os.path.join(out_dir, "log.txt"), "w")
  msg = ["Fetching GSE %s " % gse_id]
  if gpl_id is not None:
    msg.append("[GPL %s] " % gpl_id)
  timestamp = time.strftime("%a, %d %b %Y %H:%M:%S")
  script = os.getcwd()
  msg.append("on %s using %s Version %s" % (timestamp, script, VERSION))
  report("".join(msg), fp_log)
  report("Using EQTLFilter only, default parameters", fp_log)

  # Create GSE object.
  g = GSE(gse_id, platform_id=gpl_id)

  # If g is a super study, fetch all sub studies
  if g.type == "SUPER":
    report("%s is a super study with %d children." % (g, len(g.substudies)), fp_log)
    for gsub in g.substudies.values():
      # Skip substudies that are not eQTL
      if gsub.type != "eQTL":
        report("%s is type %s. Skipping..." % (gsub, gsub.type), fp_log)
        continue
      write_study(gsub, fp_log, out_dir)
  # Otherwise, simply fetch G itself.
  else:
    report("%s is a child study. Fetching it directly..." % (g), fp_log)
    write_study(g, fp_log, out_dir)


def write_study(gse, fp_log, out_dir=""):
  """Write a filtered GSE matrix to a new file.

  Args:
    gse: geo.GSE non-super study instance
    fp_log: [*str] open writable file pointer for logging
    out_dir: str of output directory
  """
  filename = "%s.%s.%s.tab" % (gse.id, gse.platform.id, gse.type)
  fp = open(os.path.join(out_dir, filename), "w")
  report("Writing %s to file %s with default EQTLFilter..." % (gse, filename), fp_log)
  
  filt2 = EQTLFilter(gse)
  n_lines = 0
  for row in filt2.get_rows():
    n_lines += 1
    # DO NOT Skip headers.
    # If this is the first line, print a "#" to indicate that this is a header line
    if n_lines == 1:
      fp.write('#')
    
    # First 5 columns: 'ID_REF', 'GENE_SYMBOL', 'NUM_VALUES', 'MEAN', 'STD'
    #   remove all but 'GENE_SYMBOL' column which should be made first column
    row = [row[1]] + row[5:]
    # replace all "None" with empty string
    for i, v in enumerate(row):
      if v == "None":
        row[i] = ""
    fp.write("\t".join(row) + "\n")
  fp.close()
    
    
if __name__ == "__main__":

  try:
    gse_id = sys.argv[1]
    options = dict(map(lambda s: s.split('='), sys.argv[2:]))
  except:
    print USE_MSG
    raise
  main(gse_id, **options)

