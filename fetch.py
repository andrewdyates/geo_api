#!/usr/bin/python
"""Write a single, filtered GSE study in MINE format to STDOUT. 
Script wrapper for geo_api. 
Does not handle pseudo-super or super studies.
"""
USE_MSG = "USE: fetch.py GSE_ID [GPL_ID]"

from geo import *
from filter import *
import sys


def fetch(gse, gpl=None):
  g = GSE(gse, platform_id=gpl) 
  filt2 = EQTLFilter(g)

  n_lines = 0
  for row in filt2.get_rows():
    n_lines += 1
    # Skip headers
    if n_lines == 1: 
      continue
    # First 5 columns: 'ID_REF', 'GENE_SYMBOL', 'NUM_VALUES', 'MEAN', 'STD'
    #   remove all but 'GENE_SYMBOL' column which should be made first column
    row = [row[1]] + row[5:]
    # replace all "None" with empty string
    for i, v in enumerate(row):
      if v == "None":
        row[i] = ""
    print "\t".join(row)
    
    
if __name__ == "__main__":
  if len(sys.argv) == 2:
    # only GSE id
    fetch(sys.argv[1])
  elif len(sys.argv) == 3:
    # GSE id, plus GPL platform for psuedo-studies
    fetch(sys.argv[1], sys.argv[2])
  else:
    print USE_MSG
    sys.exit(1)
