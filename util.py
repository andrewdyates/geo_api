#!/usr/bin/python
from geo import *

def get_gpl(gpl_brief=None, gpl_data=None, is_tab=True):
  """Return geo.Local_GPL object from GPL file names.

  Args:
    gpl_brief: str of file path to GPL file meta data brief (as from geo_downloader)
    gpl_data: str of file path to GPL row descriptions IN DATA ROW ORDER (as from geo_downloader)
      may either be .tab for SOFT text file format
    is_tab: bool if gpl_data is in .tab format
  Returns:
    geo_api.GPL of loaded row descriptions
  """
  gpl = LocalGPL(fname_brief=gpl_brief, fname_data=gpl_data, data_is_tab=is_tab)
  gpl.load()
  return gpl

def assert_row_alignment(n, varlist, gpl):
  """Assert that GPL row descriptions and row id list aligns to data rows."""
  assert n == len(gpl.row_desc) == len(gpl.probe_list) == len(gpl.probe_idx_map), \
      " != ".join((n, len(gpl.row_desc), len(gpl.probe_list), len(gpl.probe_idx_map)))
  for i in xrange(n):
    s = gpl.probe_list[i]
    assert s == varlist[i]
    assert gpl.probe_idx_map[s] == i
    assert s in gpl.row_desc
