#!/usr/bin/python
"""Loose test scripts for GEO Fetch."""
from geo import *
from filter import *
import sys
import subprocess


def print_test1(gse):
  g = GSE(gse)
  print g
  for name in dir(g):
    if name[0] == "_":
      continue
    print name, ":", getattr(g, name)
  print "========="
  print
  for substudy in g.substudies:
    print "----------"
    for name in dir(substudy):
      if name[0] == "_":
        continue
      print name, ":", getattr(substudy, name)


def test2(gse):
  """Dump GSE data rows to STDOUT"""
  g = GSE(gse)
  h = g.substudies["GSE25935"]
  print "getting rows for ^%s..." % h.id
  i = 0
  print h.platform.populated
  print len(h.platform.id_desc)
  # this may not be working correctly
  for row in h.get_rows():
    if i > 40:
      break
    print "\t".join(row)
    i += 1
  i = 0
  for k,v in h.platform.id_desc.items():
    if i > 10:
      break
    print k, v
    i += 1

def test4(gse):
  """Dump GSE substudy data rows to STDOUT test 2"""
  g = GSE(gse)
  i = 0
  for row in g.get_rows():
    if i > 40:
      break
    print "\t".join(row)
    i += 1
    

def test3(gse):
  """Filter the entire GSE study.

  How to customize Filters?
  """
  # see parser work for GSE19177
  # will need to upload blast to get genes for SNPs
  # THIS DOES NOT WORK
  g = GSE(gse)
  if g.type not in RECOGNIZED_STUDY_TYPES:
    raise Exception, "Unrecognized type %s." % g.type
  filt = EQTLFilter(g)
  i = 0
  for line in filt.get_rows():
    print line
    i += 1
    if i > 5:
      break
  print g.platform.row_desc["ILMN_1659893"]



def test5(gse):
  """Filter the entire GSE study.

  How to customize Filters?
  """
  print "TEST 5"
  g = GSE(gse)

  if g.type not in RECOGNIZED_STUDY_TYPES:
    raise Exception, "Unrecognized type %s." % g.type
  
  i = 0
  for row in g.get_rows():
    print row
    i += 1
    if i > 5:
      break
  print "what's going on?"
    

def test_plot():
  """Filter the entire GSE study.

  How to customize Filters?
  """
  g2 = GSE("GSE28893")
  filt2 = EQTLFilter(g2)
  for line in filt2.get_rows():
    pass
  
  g = GSE("GSE25935")
  filt = EQTLFilter(g)
  for line in filt.get_rows():
    pass

  x = range(len(filt.row_stats))
  y1 = [d['std'] for d in filt.row_stats.values()]
  x2 = range(len(filt2.row_stats))
  scale=len(x)/float(len(x2))
  xx2 = map(lambda x: x*scale, x2)
  
  y2 = [d['std'] for d in filt2.row_stats.values()]
  y1.sort()
  y2.sort()

  pylab.clf()
  pylab.cla()
  pylab.plot(x, y1, 'b', label="GSE28893")
  pylab.plot(xx2, y2, 'r', label="GSE28893")
  pylab.show()
  
def fetch(gse="GSE25935"):
  g = GSE(gse) 
  filt2 = EQTLFilter(g)
  i =0
  for row in filt2.get_rows():
    # yield column headers?
    i += 1
    print "\t".join(row)
    if i >= 10:
      break
  
        
    
if __name__ == "__main__":
  if len(sys.argv) >2:
    fetch(sys.argv[1])
  else:
    fetch()
