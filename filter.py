#!/usr/bin/python
"""Filter cleaned rows of GSE series matrix data.

Environment variables:
TMP_DIR: path to temporary directory

TODO: Configure filter parameters
"""
import os
import pylab
import random
import time

import geo
from logger import Log

class MalformedFilterError(Exception):
  pass

class Filter(object):
  pass

def temp_file_name(name=""):
  """Return string of a full path to a new temporary file name.

  Args:
    name: str of name to include in temporary file name.
  Returns:
    str of full path to new file name in TMP_DIR
  """
  t = int(time.time())
  x = random.randint(0, 1000000)
  filename = "tmp.%s.%d.%d" % (name, t, x)
  if 'TMP_DIR' in os.environ:
    path = os.environ['TMP_DIR']
  else:
    # The temporary, working directory should also be environment sensitive.
    raise Exception, "Set environment variable TMP_DIR to some path."
  return os.path.join(path, filename)

def get_float(s):
  """Return string as float or as None if it cannot be converted.

  Args:
    s: str, presumably decimal float or else some representation of "None"
  Returns:
    float or None
  """
  try:
    x = float(s)
  except ValueError:
    return None
  else:
    return x

def merge_titles(values):
  """Return merged values as a ';' delimited concatenation.

   Args:
    values: [str] of column values in sample equivalence class
  Returns:
    str: `values` sorted <= and joined together by the ';' character
  """
  return ';'.join(sorted(values))

def merge_floats(values):
  """Return merged values as an average of floats (or None).

  Args:
    values: [str] of column values in sample equivalence class
  Returns:
    float: average of numeric values or None if no numeric values
  """
  v_sum = 0
  n_sum = 0
  for v in values:
    # Try to convert all column values to floats.
    try:
      x = float(v)
    # Ignore values that cannot be converted to a float.
    except ValueError:
      pass
    # Increment aggregate value
    else:
      v_sum += x
      n_sum += 1
      
  # Compute and return average value if it exists, else return None
  if n_sum <= 0:
    return None
  else:
    s = v_sum / n_sum
  return s


class EQTLFilter(Filter):
  """Filter a compiled GSE dataset.
   Consider adding a "read_data" function; I may not need to emit lines.

  Attributes:
    gse: geo.GSE populated study data instance 
    col_titles: [str] of column titles of filtered data matrix
    col_map: {int:set(int)} of disjoint column equivalence classes
    rows_filtered: [str] of line ids filtered for missing a gene symbol
    rows_per_gene: {str: set(str)} of gene symbols to row ids
    row_stats: {str: {'num_values': int, 'mean': float, 'std': float}} per rowID
      note: only rows not in rows_filtered are in row_stats
  """
  def __init__(self, gse, merge_cols=True, percentile=.75):
    """Initialize filter. Requires populated gse.

    Args:
      gse: GSE instance associated with row_iter
      merge_cols: bool if to merge columns if able
      percentile: float 0<x<=1 of top percent by std to keep
    """
    # 1. Require that GSE is populated and is of correct type.
    # ==========
    if not gse.populated:
      raise geo.NotPopulatedError, "%s must be populated to filter rows." % gse
    if gse.type != "eQTL":
      raise geo.StudyTypeMismatch, "%s must be type 'eQTL', not '%s'." % \
        (gse, gse.type)

    # 2. Set Attributes.
    # ==========
    self.gse = gse
    self.col_titles = self.gse.col_titles[:]
    self.col_map = None
    self.rows_filtered = []
    self.rows_per_gene = {}
    self.row_stats = {}
    self.merge_cols = merge_cols
    self.percentile = percentile
    
    # 3. Get column map for column merging.
    # ==========
    n_samples = len(self.gse.samples)
    n_uniques = len(self.gse.subject_gsms)

    # If there are more samples than unique subjects, then create column map.
    if self.merge_cols and n_samples > n_uniques:
      self.col_map = self._make_col_map()
      rx_str = self.gse.parameters['rx_gsm_subject_str']
      Log.info(("Created column merge map for %s (%d samples to %d subjects)" +\
        " with rx '%s'") % \
        (self.gse, n_samples, n_uniques, rx_str))
      # Verify that column merge map is reasonable (num uniques + 1 for ID column)
      if len(self.col_map) != n_uniques + 1:
        Log.warning("Column merge map has %d classes, expected %d in %s." % \
                    (len(self.col_map), n_uniques, self))
        
    # No column merging scheme can exist. Do not create a col_map.
    else:
      # Retrieve the regular expression used
      rx_str = self.gse.parameters['rx_gsm_subject_str']
      Log.info("No column merge map created for %s using rx '%s'. Merge_cols flag is %s" % \
        (self.gse, rx_str, self.merge_cols))

  def __repr__(self):
    return "[EQTLFilter => %s (%d)]" % (self.gse, id(self))

  def get_rows(self):
    """Return filtered row iterator.
    CLEAN THIS UP
    It may be best to break this into multiple filters?
    Fix to return [str]
    
    Returns:
      *[str] of filtered rows of data split by columns
    """
    Log.info("Initiated filter %s for rows of %s" % (self, self.gse))
    if self.col_map:
      Log.info("self.col_map exists. Merge %d to %d columns for %s" % \
               (len(self.col_titles), len(self.col_map), self))
    else:
      Log.info("No col_map. Will not merge %d columns for %s." % \
               (len(self.col_titles), self))

    # 0. Determine best gene name column in case GENE_SYMBOL does not exist.
    # ==========
    gene_symbol_name = None
    # Traverse column names in preferred order.
    for name in geo.GPL.EQTL_GENE_NAME_LIST:
      # Skip columns without assignments. Continue
      if self.gse.platform.special_cols[name] is None:
        continue
      # Choose the first column that has an acceptable assignment. Break.
      else:
        actual_column_name = self.gse.platform.special_cols[name]
        gene_symbol_name = name
        break
    # Verify that a column was chosen to identify the row.
    if gene_symbol_name:
      Log.info("Selected column '%s=>%s' to best represent gene name for %s." %\
        (gene_symbol_name, actual_column_name, self.gse.platform))
    else:
      raise MalformedFilterError, "Cannot select gene symbol column from %s" % \
        (self.gse.platform)
    
    # 1. Update column titles accounting for merged columns.
    # ==========
    if self.col_map:
      self.col_titles = self._merge_cols(self.col_titles, merge_titles)
      
    # Insert generated column titles (AFTER merging columns)
    # self.col_titles[0] should always be "ID_REF"
    col_titles_prefix = ["ID_REF", gene_symbol_name, "NUM_VALUES", "MEAN", "STD"]
    self.col_titles = col_titles_prefix + self.col_titles[1:]
    Log.info("Added %s, NUM_VALUES, MEAN, STD to col titles for %s." %\
             (gene_symbol_name, self))
             
    # Open new temporary file. XXX RENAME
    filepath = temp_file_name("%s.rowmerge" % self.gse.id)
    fp_out = open(filepath, "w")

    # 2: @DATAPASS 1: Merge columns, add gene symbol, filter non-genes.
    # ==========
    Log.info(("Started filter 1 in %s for %s: find and add gene, merge cols. " +
             "(This may take a while.)") % (self, self.gse))
      
    num_rows = 0
    for row in self.gse.get_rows():
      # TODO: Add status reporting to console
      num_rows += 1

      # Determine gene symbol for this row. Filter if no gene symbol exists.
      row_id = row[0] # Row ID should always be the first entry in a row.
      gene_sym = self.gse.platform.get_column(row_id, gene_symbol_name)
      if not gene_sym:
        self.rows_filtered.append(row_id)
        continue # skip this row
      else:
        self.rows_per_gene.setdefault(gene_sym, set()).add(row_id)
      
      # Merge columns using column mapping of series matrix columns.
      # Also, transform row into "floats" and None
      if self.col_map:
        # XXX_merge_cols is slow, perhaps due to float conversions.
        row = self._merge_cols(row, merge_floats)
      else:
        row = map(get_float, row)

      # Compute mean and standard deviation of all non-ID columns
      # check for None specifically since a valid value could be 0
      filtered_row = filter(lambda x: x is not None, row[1:])
      std = pylab.std(filtered_row)
      mean = pylab.mean(filtered_row)
      num_values = len(filtered_row)
      # Store row statistics
      self.row_stats[row_id] = \
        {'num_values': num_values, 'mean': mean, 'std': std}

      # Insert (gene_sym, size, mean, std) into second column
      row = [row_id , gene_sym, num_values, mean, std] + row[1:]

      # Write row to temporary file.
      # TODO: I may want to compress my row by converting it to a pickle.
      # pickling a list of floats uses 2/3 space and takes 1/2 compute time.
      fp_out.write("\t".join(map(str, row)))
      fp_out.write("\n")
    fp_out.close()

    # Log results of filter pass 1
    # ==========
    n = len(self.rows_filtered)
    n_gene_rows = num_rows-n
    mean_rows_per_gene = float(num_rows-n)/len(self.rows_per_gene)
    
    if num_rows != self.gse.est_num_row:
      Log.warning("Num rows read(%d) not num rows expected(%d) for %s" % \
                  (num_rows, self.gse.est_num_row, self))
    Log.info(("Filter 1 complete for %s. " + \
      "%d of %d (%.2f%%) rows removed for no gene symbol. %d rows remain.") % \
      (self, n, num_rows, (n/float(num_rows))*100, n_gene_rows))
    Log.info("Number of unique genes: %d, %.1f mean num rows per gene." % \
      (len(self.rows_per_gene), mean_rows_per_gene))

    # 3: Choose representative genes from self.row_stats and self.rows_per_gene
    # ==========
    # select all rows for a gene. If a gene 
    selected_row_ids = []
    for gene, row_ids in self.rows_per_gene.items():
      # If only a single row for this gene exists, choose it.
      if len(row_ids) == 1:
        best_row_id = row_ids.pop()
      # Else, choose row with the highest mean value. 
      else:
        s = sorted(row_ids, key=lambda x: self.row_stats[x]['mean'])
        best_row_id = s[-1]
      # Add this row_id to the accepted list
      selected_row_ids.append(best_row_id)

    n_single_gene_rows = len(selected_row_ids)
    Log.info("Selected %d of %d rows for %d genes by maximum row mean." % \
      (n_single_gene_rows, n_gene_rows, len(self.rows_per_gene)))

    # Sort row_ids by row standard deviation in decreasing order.
    selected_row_ids.sort(key=lambda x: self.row_stats[x]['std'], reverse=True)
    
    # Select top percentile by std. Convert type to set for easier membership tests.
    x = int(len(selected_row_ids)*self.percentile)
    selected_row_ids = set(selected_row_ids[:x])
    threshold_num_rows = len(selected_row_ids)
    assert(x == threshold_num_rows)
    Log.info("Selected top %d%% of rows (%d of %d) by standard deviation." % 
      (self.percentile*100, threshold_num_rows, n_single_gene_rows))
      
    # FINAL PASS: YIELD FILTERED LINES
    # ===========
    # Open temporary file generated in first pass.
    fp = open(filepath, "r")

    # Yield (modified) column titles.
    yield self.col_titles[:]
    
    # For each line, only yield if the row_id is in the selected_row_ids list.
    num_yielded_rows = 0
    for line in fp:
      row = line.strip().split("\t")
      row_id = row[0]
      if row_id in selected_row_ids:
        num_yielded_rows += 1
        yield row

    # All lines yielded. Check number of lines yielded with expected value.
    if num_yielded_rows != threshold_num_rows:
      Log.warning("%d yielded rows != %d expected number of rows." % \
        (num_yielded_rows, threshold_num_rows))
    else:
      Log.info("Filter complete. yielded %d rows." % (num_yielded_rows))
      
      
  def _merge_cols(self, row, f_merge):
    """Return column-merged row.

    Args:
      row: [value] of column values
      f_merge: function([values]) with which to merge column class values
    Return:
      [value] of merged column values, len <= len(`row`)
    """
    # Verify that the row aligns with the column map and that the map exists.
    if not self.col_map:
      raise MalformedFilterError, \
        "_merge_cols() called, but col_map does not exist for %s" % self
    
    # Merge columns per row.
    new_row = []
    consumed_cols = set()
    # For each column in this row from the left:
    for i, col in enumerate(row):
        
      # Ignore consumed columns
      if i in consumed_cols:
        continue
        
      # Trivially add columns with an undefined equivalence class.
      if self.col_map[i] is None:
        consumed_cols.add(i)
        new_row.append(col)
        continue

      # Merge all values of equivalence class and append aggregate value.
      s = [row[j] for j in self.col_map[i]]
      v = f_merge(s)
      new_row.append(v)
      # Consume all columns in equivalence class.
      consumed_cols.update(self.col_map[i])

    # Verify that all columns in row have been consumed
    if len(row) != len(consumed_cols) or len(new_row) != len(self.col_map):
      raise MalformedFilterError, \
        ("Mismatch in column merge for %s. " % self) + \
        "#cols in: row=%d, consumed=%d, new_row=%d, col_map=%d" % \
        (len(row), len(consumed_cols), len(new_row), len(self.col_map))
          
    return new_row
    
  def _make_col_map(self):
    """Return column equivalence class map used in per-row column merging.

    Returns:
      {int, set(int)} of (unmodified) row id equivalence classes
    """
    # Set of columns already mapped.
    consumed_cols = set()
    # int=>set(int) of first col num representative to its equivalence class.
    col_map = {}
    for i, title in enumerate(self.col_titles):
      # Column is not a sample. Its equivalence class is undefined.
      if title not in self.gse.samples:
        col_map[i] = None
        consumed_cols.add(i)
      else:
        # Select sample's subject from GSM object given this column title.
        subject = self.gse.samples[title].subject
        # For each GSM ID mapped to this subject:
        for gsm in self.gse.subject_gsms[subject]:
          # Get corresponding col num j for this gsm. j may equal i.
          j = self.col_titles.index(gsm)
          # Ignore consumed columns.
          if j in consumed_cols:
            continue
          # Map col j to col i. Add j to the list of consumed columns.
          col_map.setdefault(i, set()).add(j)
          consumed_cols.add(j)

    return col_map
