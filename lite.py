#!/usr/bin/python
"""Simplistic wrappers for locally loaded GEO objects.
"""
import re
import collections

RX_PLATFORM = re.compile("^\^PLATFORM = (\w+)")
RX_HEADER = re.compile("^#([^=]+?) = ?(.*)")
RX_ATTR = re.compile("^!Platform_([^=]+?) = ?(.*)")
HEAD_END_LINE = "!platform_table_begin"
TABLE_END_LINE = "!platform_table_end"

class GPL_Lite(object):

  def __init__(self, filename):
    self.cols = collections.OrderedDict()
    self.rows = collections.OrderedDict()
    
    fp = open(filename)
    # GPL ID
    self.id = RX_PLATFORM.match(fp.next()).group(1)
    # Column Attribute Descriptions
    for line in fp:
      m = RX_HEADER.match(line)
      if not m: break
      key, value = m.groups()
      self.cols[key] = value
    assert line[:-1] == HEAD_END_LINE
    # Column titles
    assert self.cols.keys() == fp.next().strip('\n').split('\t')
    # Data rows
    for line in fp:
      line = line.strip('\n')
      if line == "": continue
      if line == TABLE_END_LINE: break
      row = line.split('\t')
      self.rows[row[0]] = dict(zip(self.cols[1:], row[1:]))
        
    
