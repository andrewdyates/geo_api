#!/usr/bin/python
"""Classes for filtering GSE datasets."""

class RowFilter(object):
  def __init__(self, gse):
    """Initialize filter.

    Args:
      gse: geo.GSE object associated with this filter.
    """
    self.gse = gse
    
  def get(self, row):
    return row
  def _insert(self, row, value):
    row.insert(self.COL_NUM, value)

class RemoveID(RowFilter):
  def __init__(self, gse):
    """Initialize filter.

    Args:
      gse: geo.GSE object associated with this filter.
    """
    self.gse = gse
  def get(self, row):
    return row
  

class RowRun(object):
  def __init__(self, gse, f_reduce, f_key=None):
    self.gse = gse
    self.f_reduce = f_reduce
    if f_key is None:
      self.f_key = lambda r: r[0]
    self.results = {}

  def add(row):
    key = self.f_key(row)
    value = self.f_reduce(row)
    self.results[key] = value

  def stop():
    # Do computations like 25 percentile variance threshold
    pass

  def fetch(row):
    return self.results[f_key(row)]


class InsertGene(RowFilter):
  """Add gene symbol to row."""

  # All rows should be begin with its row id.
  COL_NUM = 0
  def __init__(self, gse):
    self.gse = gse
    # Original id column.
    self.input_id_col = 0
    
    # Verify that platform has "gene symbol" column.
    # This logic may warrant something more sophesticated in the future.
    # For example, it may be provided by the GPL object.
    if "GENE_SYMBOL" not in self.platform.col_titles:
      raise ValueError, "GENE_SYM not in %s. GPL col titles are: %s" \
         % (self.gse.platform, self.gse.platform.col_titles)
    # update gse col_titles to include new gene_sym row
    self._insert(gse.col_titles, "GENE_SYM")

  def get(self, row):
    """Return modified row.

    Args:
      row: [str] of data row in expected filter order state.
    Returns:
      [str] of modified row or None if row has been removed.
    """
    row_id = row[self.input_id_col]
    desc = self.platform.id_desc[row_id]
    try:
      sym = desc["GENE_SYM"]
    except KeyError:
      sym = ""
    self._insert(row, sym)
    return row


class RemoveNonGene(RowFilter):

  def __init__(self, gse):
    self.gse = gse

    # Get gene_symb col from GSE (and verify it exists)
    if "GENE_SYMBOL" not in self.col_titles:
      raise ValueError, "GENE_SYM not in %s (maybe this filter is misordered?) GSE col titles are: %s" \
         % (self.gse, self.gse.col_titles)
    self.input_sym_col = self.gse.col_titles.index("GENE_SYMBOL")
    
  def get(self, row):
    if not row[self.input_sym_col]:
      return None
    else:
      return row
    
class MergeCols(RowFilter):
  pass
