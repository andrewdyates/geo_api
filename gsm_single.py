#!/usr/bin/python
"""Parse a GSM "full" SOFT text file downloaded from GEO.

Includes sample data.
"""
import re

RX_TITLE = re.compile("^\^SAMPLE = (\w+)")
RX_SAMPLE = re.compile("^!Sample_(\w+) = (.+)$")
RX_DESC = re.compile("^#([^=]+?) = ?(.*)")
S_BEGIN = "!sample_table_begin"
S_END = "!sample_table_end"

# Create download URL to SOFT text "full" using URL_PTN % GSM_ID
URL_PTN = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=full"


class GSM_Lite(object):
  def __init__(self, fp):
    self.id = None
    self.attr = {}
    self.col_desc = {}
    self.row_ids = []
    self.rows = {}

    # title
    self.id = RX_TITLE.match(fp.next().strip('\n\r')).group(1)
    # headers
    for line in fp:
      line = line.strip('\n\r')
      if line == S_BEGIN:
        break
      m = RX_SAMPLE.match(line)
      if m:
        self.attr.setdefault(m.group(1), []).append(m.group(2))
      else:
        m = RX_DESC.match(line)
        if m:
          self.col_desc.setdefault(m.group(1), []).append(m.group(2))
    self.split_derivatives()
    # data table column headers
    self.col_headers = fp.next().strip('\n\r').split('\t')
    for line in fp:
      line = line.rstrip('\n\r')
      if line == S_END:
        break
      row = line.split('\t')
      self.row_ids.append(row[0])
      self.rows[row[0]] = row[1:]
      

  def split_derivatives(self):
    """Attempt to split multi-attributes assigned to a single attribute key."""
    for k, v in self.attr.items():
      if len(v) > 1:
        new_attrs = {}
        old_values = set()
        for x in v:
          q = [s.strip() for s in x.split(':')]
          if len(q) == 2:
            new_attrs.setdefault("%s:%s"%(k,q[0]), set()).add(q[1])
          else:
            old_values.add(x)
        self.attr[k] = old_values
        self.attr.update(new_attrs)
      
      
