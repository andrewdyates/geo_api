#!/usr/bin/python
"""GEO (Genetic Expression Omnibus) self-contained objects.

Andrew D. Yates

TODO: 
  add support for miRNA expression type studies
  add support for SNP genomic sequencing type studies
"""
import re
import csv
import sys

from download import Download
# patched, local version of gzip from Python 3 to handle http streams, stream closes
from download import Gzipper

from logger import Log

RECOGNIZED_STUDY_TYPES = set(["eQTL", "SNP", "SUPER"])
# Default "no-merging" sample title pattern
DFT_RX_TITLE_STR = "(.*)()"
# Default GSE mutable parameters
GSE_PARAMETERS_DEFAULT = {'rx_gsm_subject_str': DFT_RX_TITLE_STR}


class MalformedDataError(Exception):
  pass
class SuperStudyAccessError(Exception):
  pass
class NotPopulatedError(Exception):
  pass
class StudyTypeMismatch(Exception):
  pass


# TODO: this should be loaded from an external settings file.
GSE_SETTINGS = {
  'GSE25935': {'rx_gsm_subject_str': "([^_]+)(?:_rep(\d+))?" },
}

# TODO: this should be loaded from an external settings file.
GPL_SETTINGS = {
#  'GPL4133': {'gene_symbol_title': "symbol" },
}

def keyword_score(title, desc, keywords):
  """Return match score for column title:desc given keyword list.
  
  Args:
    title: str of column title
    desc: str of column description
    keywords: list of keywords to match
  Returns:
    int of weighted keyword matches
  """
  score = 0
  title_words = set(re.split("[^a-z]+", title.lower()))
  desc_words = set(re.split("[^a-z]+", desc.lower()))
  score += 2.0 * len(keywords & title_words)
  score += 1.0 * len(keywords & desc_words)
  return score


class GSE(object):
  """A GEO genetic study.

  Attributes:
    id: str of GSE ID like GSE\d+
    attr: {str: [str]} of key=>values from GSE Brief, values order by appearance
    type: str of study data type in RECOGNIZED_STUDY_TYPES
    super_id: str of super study GSE ID or None if has no super study
    substudies: {str:GSE} of gse_id=>GSE substudy objects
    samples: {str: GSM} of GSM_ID=>obj fetched sample definitions
    platform: GPL of platform definition object
    parameters: {str: obj} of user-specified metadata not derived from GSE brief
      Recognized parameters:
        'rx_gsm_subject_str': str of RX to split GSM titles to (subject, repetition)
    populated: bool if contained GPLs and GSMs are populated. If super, then is
      if all substudies are populated
    pseudo: bool if this is an implicit super or substudy due to multiple GPLs
    subject_gsms: {str: [str]} of title_subject => [GSM ID]
    est_num_row: int of number of data rows expected
    col_titles: [str] of column titles in order from series matrix files
    selected_platform_id: str of specified GPL of this pseudo substudy
    _id_col_idx: [int] of columns representing ID_REF

  In GSE Super Studies (type=super), related GPL and GSM instances reside
    in their respective substudies, not in the Super Study's instance itself.
  """
  # Related to GSE Brief SOFT text files.
  PTN_GSE_BRIEF = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%(id)s&targ=self&view=brief&form=text"
  RX_SERIES = re.compile("\^SERIES = (\w+)")
  RX_HEADER = re.compile("^!Series_(\S+) = ?(.*)$")
  RX_SUBSTUDY = re.compile("^SuperSeries of: (GSE\d+)$")
  # Related to GSE series matrix SOFT text files.
  PTN_GSE_SERIES_DIR = "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/%(id)s/"
  RX_SAMPLE_HEADER = re.compile("^!Sample_(\w+)\t(.+)$")
  HEAD_END_LINE = "!series_matrix_table_begin"

  EQTL_TYPE_LINES = set([
    "Expression profiling by array",
    ])
  SNP_TYPE_LINES = set([
    "SNP genotyping by SNP array",
    "Genome variation profiling by SNP array",
    ])
  # keep newline for easier detection
  END_LINE = "!series_matrix_table_end\n"

  def __init__(self, gse_id, super_id=None, custom_parameters=None, \
               populate=True, platform_id=None):
    """Initialize GSE object.
    Warning: if "populate" is False, the __init__ blocks serial execution
      to populate itself which may be slow if from disk or from FTP.

    Args:
      gse_id: str of valid GSE id like GSE\d+
      super_id: str of valid GSE id like GSE\d+ if this is a substudy
      custom_parameters: {str:obj} of user-specified per-instance GSE parameters
      populate: bool to populate this object immediately on instantiation
      platform_id: str of one of multiple GPLs, flags "pseudo" if not None
    """
    Log.info("Initializing new GEO object (%d) for GSE_ID: %s, GPL_ID: %s" % \
      (id(self), gse_id, platform_id))
    self.id = gse_id
    self.attr = {}
    self.type = None
    self.super_id = super_id
    # i.e., has the big dataset for this GSE been read yet?
    self.populated = False
    # is this a "fake" super study or sub study due to multiple GPLs?
    self.pseudo = (platform_id is not None)
    # selected_platform_id is only used for pseudo substudies
    self.selected_platform_id = platform_id
    if self.pseudo:
      Log.info("Flagged %s (%d) as pseudo during init, platform_id = %s" %\
        (gse_id, id(self), platform_id))
    
    self.substudies = {}
    self.platform = None
    self.samples = {}
    self.subject_gsms = {}
    self.col_titles = []
    self.est_num_row = None

    # Set mutable parameters to default values
    self.parameters = GSE_PARAMETERS_DEFAULT.copy()
    # Update parameters with custom parameters from settings.
    self.parameters.update(GSE_SETTINGS.get(gse_id, {}))
    # Update parameters with custom parameters specified by user.
    if custom_parameters:
      self.parameters.update(custom_parameters)
      
    # Populate attributes with abbreviated metadata. 
    self._load_brief()
    # ....
    # Fully populate self and substudies.
    if populate:
      self.populate()

  def populate(self):
    """Load and set values to all object attributes."""
    Log.info("Populating %s..." % self)
    if self.populated:
      Log.warning("%s is already populated." % self)
      return
    
    if self.type == "SUPER":
      Log.info("Populating SuperStudy %s." % self)
      for study in self.substudies.values():
        study.populate()
    else:
      # Do the actual hard work of loading and setting attributes.
      self._populate()
    Log.info("%s Populated." % self)
      
    self.populated = True

    
  def _populate(self):
    """Populate self from series matrix headers downloaded from FTP."""
    # 2. Get list of series matrix data files
    # ==========
    ftp_files = self._get_ftp_files()

    # TODO: if pseudo, only load series matrix files with same GPL

    # 3. Open and decompress each file as a list of file pointers.
    # ==========
    Log.info("Fetching %d file(s): %s" % (len(ftp_files), ftp_files))
    fps = self._open_ftp_files(ftp_files, finalize=False)

    # 4. Populate GSM sample objects from series matrix files. 
    # ==========
    for fp in fps:
      Log.info("Populating GSMs from file pointer %s" % (fp))
      self._populate_gsms(fp)
      Log.info("GSMs Populated from file pointer, closing %s" % (fp))
      # fp now points to first line after header
      
    # Check that GSM samples have been populated.
    # Sometimes, not all samples are included in substudies.
    empty_keys = filter(lambda x: not self.samples[x].populated, self.samples)
    n_empty = len(empty_keys)
    if n_empty > 0:
      # If this is a pseudo-study, then filter samples expected from the super study
      #   but which are not actually populated from the series matrix file.
      if self.pseudo:
        for key in empty_keys:
          del self.samples[key]
        Log.info("Filtered %d of %d samples from series matrix for %s." % \
          (len(self.samples), n_empty+len(self.samples), self))
      else:
        Log.warning("Not all samples in %s (only %d of %d) have attributes." % \
          (self, n_empty, len(self.samples)))

    # 5. Populate column titles given fp pointing after header. Close fp.
    # ==========
    for fp in fps:
      self._populate_col_titles(fp)
      fp.close()

    # 6. Estimate number of data rows from first sample attribute "data_row_count"
    # ==========
    if len(self.samples) > 0:
      # get all sample's data_row_counts
      try:
        s = [int(q.attr['data_row_count'][0]) for q in self.samples.values()]
        self.est_num_row = max(s)
      except Exception, e:
        Log.warning("Could not estimate number of data rows for %s: %s" % \
          (self, e))
        

    # 7. Log population result and statistics
    # ==========
    # warn about unexpected population results
    if len(self.col_titles) == 0:
      Log.warning("0 column titles populated for %s.")
    if len(self.samples) == 0:
      Log.warning("0 Samples populated for %s.")
    if len(self.subject_gsms) == 0:
      Log.warning("0 Subjects populated for %s.")
      
    Log.info(("Populated substudy %s. Cols=%d, Samples=%d, Subjects=%d, " 
             + "Expected Rows=%d") % \
             (self, len(self.col_titles), len(self.samples), \
              len(self.subject_gsms), self.est_num_row))


  def _open_ftp_files(self, ftp_files, finalize=True):
    """Return list of open, decompressed file pointers from file list.

    Args:
      ftp_files: [FTPFile] of files to open for downloading
      finalize: bool to open files for complete (rather than partial) download
    Returns:
      [*str] of open file pointer-like objects for each file in `ftp_files`
    """
    fps = []
    for ftp_file in ftp_files:
      handle = Download(ftp_file.url, finalize=finalize)
      http_fp = handle.read()
      if ftp_file.compressed:
        # closing this file pointer should close the underlying buffer.
        zip_fp = Gzipper(fileobj=http_fp, mode='r')
        fps.append(zip_fp)
      else:
        fps.append(http_fp)
    return fps

  def _get_ftp_files(self):
    """Return list of FTP page objects from ftp page.

    Returns:
      [FTPFile] from GSE_SERIES FTP page.
    """
    ftp_files = []
    root_url = self.PTN_GSE_SERIES_DIR % {'id': self.id}
    handle = Download(root_url)
    http_fp = handle.read()
    for line in http_fp:
      f = FTPFile(root_url, line.strip())
      ftp_files.append(f)
    http_fp.close()

    # If this is pseudo-study, an FTP file may refer to a sibling study.
    # Filter all ftp files that do not include this study's GPL ID 
    if self.pseudo:
      n = len(ftp_files)
      s = self.selected_platform_id.lower()
      ftp_files = filter(lambda x: s in x.filename.lower(), ftp_files)
      Log.info("Selected %d of %d ftp files by '%s' in filename for %s." % \
         (len(ftp_files), n, s, self))
    
    return ftp_files
    
  def __repr__(self):
    if not self.pseudo:
      return "[GEO %s (%d)]" % (self.id, id(self))
    else:
      # Include platform in name to distinguish pseudo sub-studies.
      return "[GEO %s-%s (%d)]" % (self.id, self.selected_platform_id, id(self))
    
  def _load_brief(self):
    """Fetch, open, and load GSE brief metadata into attributes."""
    url = self.PTN_GSE_BRIEF % {'id': self.id}
    handle = Download(url)
    http_fp = handle.read()
    self._parse_brief(http_fp)
    http_fp.close()

  def _parse_brief(self, fp):
    """Parse GSE text brief from GEO website.

    Args:
      fp: iter=>str of "!header" format lines
    """
    # 1. Verify fp by consuming the first line for its study id.
    # ==========
    line = fp.next().strip()
    try:
      series_id = self.RX_SERIES.match(line).group(1)
    except:
      raise MalformedDataError, \
        "Cannot recognize GSE ID in line %s while fetching GEO ID '%s'" \
        % (line, self.id)
    if series_id != self.id:
      raise MalformedDataError, \
          "GSE ID %s differs from requested ID %s." % (series_id, self.id)
    Log.info("Successfully fetched brief for %s from line '%s'." % (self, line))

    # 2. Interpret next lines as "!" prefixed attributes.
    # ==========
    for line in fp:
      line = line.strip()
      if line == "":
        continue
      try:
        key, value = self.RX_HEADER.match(line).groups()
      except:
        raise MalformedDataError, "Cannot parse header line '%s' in %s" % \
          (line, self)
      key, value = key.strip(), value.strip()
      # Add header to attribute dict.
      self.attr.setdefault(key, []).append(value)

    # 3. Create substudies (if any) from header "Series_relation" value.
    # ==========
    try:
      relations = self.attr["relation"]
    except KeyError:
      # This study still may have multiple GPLs.
      pass
    else:
      for line in relations:
        substudy_m = self.RX_SUBSTUDY.match(line)
        if substudy_m:
          gse_id = substudy_m.group(1)
          # Do not immediately populate substudies.
          self.substudies[gse_id] = GSE(gse_id, self.id, populate=False)

    # 4. Identify superstudies
    # ==========
    # If any substudies exist, then this is type "SUPER"
    if self.substudies:
      self.type = "SUPER"
      Log.info("%s assigned type '%s' because it contains substudies %d %s." % \
               (self, self.type, len(self.substudies), self.substudies))
    
    # If multiple platfomms for the same study, this study is "pseudo super"
    # ==========
    if len(self.attr["platform_id"]) > 1 and self.selected_platform_id is None:
      self.type = "SUPER"
      self.pseudo = True
      n_platforms = len(self.attr["platform_id"])
      Log.info("No platform_id specified for %s, but platform is ambiguous." % self)
      Log.info("%s assigned type '%s' because it contains %d platforms %s." % \
               (self, self.type, n_platforms, self.attr["platform_id"]))
      Log.info("Set pseudo=True for %s. Creating %d pseudo-substudies..." % (self, n_platforms))
      # Create pseudo-substudies, one per GPL
      for line in self.attr["platform_id"]:
        gse_id = self.id
        gpl_id = line.strip()
        self.substudies["%s-%s" % (gse_id, gpl_id)] = \
          GSE(gse_id, self.id, populate=False, platform_id=gpl_id)

    # If this is a super study, population is complete. Exit.
    if self.type == "SUPER":
      Log.info("Population complete for superstudy %s." % self)
      return

    # ------
    # note: if superstudy, exit by now
    # ------
    # 6. Make dictionary of unpopulated GSM samples.
    # ==========
    rx_gsm_subject_str = \
      self.parameters.get("rx_gsm_subject_str", None)
    # XXX: this assumption may not hold for pseudo studies
    for gsm_id in self.attr["sample_id"]:
      self.samples[gsm_id] = GSM(gsm_id, rx_title_str=rx_gsm_subject_str)


    # 5. Create platform and determine (sub)study type
    # ==========
    type_set = set(self.attr.setdefault("type", []))

    # 5.0. Confirm that declared study types seem reasonable
    # ----------
    # Determine that there exists some overlap between declared
    #   study types and supported study types
    if not (type_set & (self.EQTL_TYPE_LINES | self.SNP_TYPE_LINES)):
      Log.warning("No type recognized in declared types for %s. Types: %s" % \
                  (self, type_set))
    # Report multiple study types in a leaf study.
    if len(type_set) > 1:
      Log.info("Multiple study type declarations %s for %s." % \
                  (self.attr["type"], self))
    elif len(type_set) == 0:                  
      Log.warning("Missing study type declaration for %s." % self)

    
    # 5.1. If this is a pseudo substudy, don't rely on the given study types from the
    #   GSE metadata.
    # ----------
    if self.pseudo:
      # Create study platform with declaring a type
      gpl_id = self.selected_platform_id
      Log.info("Attempting to determine pseudo substudy %s type from %s" % \
        (self, gpl_id))
      self.platform = GPL(gpl_id, study_type=None) # the GPL type is as yet unknown
      # Use the study type guessed by the study platform based on its definition
      self.type = self.platform.type
      # Report final study type determined.
      Log.info("%s assigned type %s from platform %s." % \
               (self, self.type, self.platform))

    # 5.2 This is not a pseudo substudy; determine type from study description
    # ----------
    else:
      # eQTL expression type study
      if self.EQTL_TYPE_LINES & type_set:
        self.type = "eQTL"
      # SNP type studies
      elif self.SNP_TYPE_LINES & type_set:
        self.type = "SNP"
      # Report unrecognized type.
      else:
        Log.warning("Unrecognized study type descriptions %s for %s." % \
          (self.attr["type"], self))
        self.type = "OTHER"
      # Report final study type determined.
      Log.info("%s assigned type %s given type description '%s'." % \
               (self, self.type, self.attr["type"]))
      # Create study platform given determined study type
      gpl_id = self.attr["platform_id"][0]
      self.platform = GPL(gpl_id, study_type=self.type) # GPL type known.

    # Populate complete for substudy.
    # ==========
    Log.info("Population complete for substudy %s." % self)

      
  def get_rows(self):
    """Yield rows of unfiltered, compiled series matrix data.
    Populate all child objects (platforms, samples) if not already populated.

    Yields:
      [str] of columns of data per row
    """
    Log.info("Yielding data rows for %s..." % self)
    
    # If this is a super series, raise an exception. 
    if self.type == "SUPER":
      raise SuperStudyAccessError, \
        "%s is a super study. Get data rows from substudies %s." % \
        (self.id, self.substudies)
    # If this is not yet populated, warn in logs and populate self.
    if not self.populated:
      Log.warning("%s get_rows() called before populating; populating." % self)
      self.populate()

    # 0. Load GPL platform row definition (if not already loaded)
    # ==========
    if not self.platform.loaded:
      self.platform.load()
    else:
      Log.warning("%s of %s already loaded." % (self.platform, self))

    # 1. Get list of series matrix data files.
    # ==========
    ftp_files = self._get_ftp_files()

    # 2. Open and decompress each file as a list of file pointers.
    # ==========
    Log.info("Fetching %d file(s): %s" % (len(ftp_files), ftp_files))
    fps = self._open_ftp_files(ftp_files)
    
    # 3. Consume GSE Series Matrix headers.
    # ==========
    for fp in fps:
      self._consume_header(fp)

    # 4. Consume GSE Series Matrix column title lines
    # ==========
    for fp in fps:
      line = fp.next()
      # ID_REF should be in the column titles. Warn if it is not.
      if "ID_REF" not in line:
        Log.warning("'%s' may not be column title line as expected for %s." % \
                    line, self)

    # 5. Read study data in parallel. Call row hook function for each line.
    # ==========
    n_rows = 0
    for row in self._yield_rows(fps):
      n_rows += 1
      yield row

    # 6. Report row yield.
    # ==========
    Log.info("Yielded %d rows of data, %d expected for %s." % \
             (n_rows, self.est_num_row, self))
    if n_rows != self.est_num_row:
      Log.warning("Yielded num rows(%d) not expected num(%d) for %s" % \
                  (n_rows, self.est_num_row, self))


  def _yield_rows(self, fps):
    """Yield a row of data from a list of parallel file iterators.

    Args:
      fps: [iter=>str] of parallel file pointers
    Yields:
      [str] of columns of values
    """
    # Generator loop: read from each fp and yield one row per iteration.
    while True:
      try:
        lines = [fp.next() for fp in fps]
      except StopIteration:
        # Verify that all fp have reached end of file.
        if not self._verify_all_fp_at_eof(fps):
          raise MalformedDataError, \
            "EOF mismatch while reading %d series matrix files for %s." % \
            (len(fps), self)
        else:
          Log.info("All %d files read to EOF for %s." % (len(fps), self))
        # OK: all file pointers stopped simultaneously. Break generator loop.
        break
        
      # check for !series_matrix_table_end end line, do not yield this line
      if "!series_matrix_table_end\n" in lines:
        continue
      
      # Merge lines into a single row
      row = self._merge_csv_row_lines(lines)
      # Warn if number of columns does not match column titles
      if len(row) != len(self.col_titles):
        Log.warning(("Parsed row of %d columns != expected %d columns. " + \
                    "GSE object: %s, lines: %s") % \
                    (len(row), len(self.col_titles), self, lines))
      # Finally, yield one combined row of data in the Generator loop
      yield row

  def _merge_csv_row_lines(self, lines):
    """Return row from list of csv lines of text. Verify that row ids align.

    Args:
      lines: [str] of adjacent csv series matrix data rows
    Returns:
      [str] row of column values from compiling `lines` in order
    """
    # Given a list of lines, split them as csv, combine, and return row.
    row_id = None
    row = []
    for line in lines:
      s = csv.reader([line], delimiter="\t").next()

      # Verify that row IDs match. Discard superfluous row IDs.
      if row_id is None:
        row_id = s[0]
        row.extend(s)
      else:
        alt_id = s[0]
        if alt_id != row_id:
          raise MalformedDataError, \
            "Row ID Mismatch: %s != %s in GSE %s" % (alt_id, row_id, self.id)
        # Do not add the redundant matching row id.
        row.extend(s[1:])
        
    return row

  @staticmethod
  def _verify_all_fp_at_eof(fps):
    """Return True if all fp throw StopIteration, i.e., all fps @ EOF

    Args:
      fps: [*iter], expected to be list of file pointers
    Returns:
      bool: if all "file pointers" in fps at EOF 
    """
    for fp in fps:
      try:
        fp.next()
      except StopIteration:
        pass # OK, EOF
      else:
        return False # FAIL: StopIteration should have been raised.
    return True # Success: no failures
    
  def _populate_col_titles(self, fp):
    """Populate column titles. 
    
    Assume that fp order is preserved, that is, the order in which column titles
      are populated now will be the order in which they will be read from later.
    
    Args:
      fp: iter=>str of file pointer pointing to column titles
    """
    line = fp.next()
    row = csv.reader([line], delimiter="\t", quotechar='"').next()
      
    # Verify that first column title is ID_REF per file.
    if row[0] != "ID_REF":
      raise MalformedDataError, \
        "ID_REF not first column title in %s for %s." % (line, self)

    # Only add the first ID_REF column title.
    if len(self.col_titles) == 0:
      self.col_titles.extend(row)
    else:
      self.col_titles.extend(row[1:])

  def _consume_header(self, fp):
    """Consume series matrix header."""
    num_lines_consumed = 0
    for line in fp:
      num_lines_consumed += 1
      if line.strip() == self.HEAD_END_LINE:
        break
    Log.info("Consumed %d lines from %s for %s" % \
             (num_lines_consumed, self, fp))

  def _populate_gsms(self, fp):
    """Populate GSM list from series matrix SOFT "!" headers.
    Modifies `fp` to point to first line after headers.

    Args:
      fp: iter of str lines of GSE "!" headers; fp points to top of file.
    """
    # sample_attrs[key] = [[str]]
    sample_attrs = {}
    
    for line in fp:
      line = line.strip()
      
      # Exit loop leaving fp pointing at first line after header.
      if line == self.HEAD_END_LINE:
        break

      # Ignore !Sample headers. Collect !Sample headers.
      m = self.RX_SAMPLE_HEADER.match(line)
      if m:
        key, values = m.groups()
        # Split line by tabs into columns
        row = csv.reader([values], delimiter="\t", quotechar='"').next()

        # Add this row to list of values for this attribute.
        sample_attrs.setdefault(key, []).append(row)

    # All header lines have been consumed.
    # Map column entries per row to GSE samples instances by GSE ID order.
    # We assume that there exists only one row for "geo_accession"
    sample_list = sample_attrs["geo_accession"][0]
    Log.info("Populating %d GSM samples for %s" % (len(sample_list), self))

    # Use position "i" to map this sample with its column of attributes.
    for i, gsm_id in enumerate(sample_list):

      # Populate sample from its corresponding row in the series matrix file.
      sample = self.samples[gsm_id]
      for key, rows in sample_attrs.items():
        for row in rows:
          sample.add_pair(key, row[i])
        sample.populated = True
            
      # Map this samples' subject to its GSE ID
      self.subject_gsms.setdefault(sample.subject, []).append(gsm_id)
      
    # Correct multiple sample attributes assigned to the same name
    for gsm_id, sample in self.samples.items():
      sample.split_derivatives()
    

class FTPFile(object):
  """An abstract data file object.

  Attributes:
    url: str of full path to file
    filename: str of file name
    compressed: bool if this file seemed to be compressed
    size: int of reported file size in bytes

  FTP line format:
  -r--r--r--   1 ftp      anonymous 52831430 Dec 28 07:34 GSE25935_series_matrix-1.txt.gz
  0:permission, 1:?, 2:protocol, 3:user, 4:size, 5:month, 6:date, 7:time, 8:name
  """
  def __init__(self, dir_url, ftp_line):
    """
    Args:
      dir_url: str of full ftp path to this file's container
      ftp_line: str of whitespace delimited remote file attributes
    """
    s = re.split("\s+", ftp_line.strip())
    self.filename = s[8]
    self.size = int(s[4])
    self.compressed = (self.filename.split(".")[-1].lower() == "gz")
    self.url = "%s/%s" % (dir_url.rstrip('/'), self.filename)

  def __repr__(self):
    return "[FTPFile %s (%d)]" % (self.filename, id(self))

  
class GPL(object):
  """A GEO Platform definition.

  Attributes:
    id: str of GEO id like GPL4133
    row_desc: {str: {str: str}} of row ids to dict of col_title=>value per row
    type: str in RECOGNIZED_STUDY_TYPES of platform type
    populated: bool if this GPL has metadata
    loaded: bool if this GPL has row descriptions
    col_titles: [str] of column titles in order they appear in GPL text matrix
    col_desc: {str:str} of column names to descriptions
    attrs: {str: [str]} of GPL attributes
    special_cols: {str:str} of special column title to its actual column title
    case_insensitive_row_id: {str:str} of row_id in lowercase to actual case
    parameters: {str: obj} of supplied metadata not derived from GPL definition
      Recognized parameters:
        all KEYWORDS column names, for example:
        'GENE_SYMBOL': str of column title containing a row's gene symbol
  """
  PTN_GPL = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%(id)s&targ=gpl&view=data&form=text"
  PTN_GPL_QUICK = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%(id)s&targ=self&view=quick&form=text"
  RX_PLATFORM = re.compile("^\^PLATFORM = (\w+)")
  # matches key, value
  RX_HEADER = re.compile("^#([^=]+?) = ?(.*)")
  RX_ATTR = re.compile("^!Platform_([^=]+?) = ?(.*)")
  HEAD_END_LINE = "!platform_table_begin"
  TABLE_END_LINE = "!platform_table_end"

  # Use these keywords to find special columns and determine study type
  # {type: {col_name: [key_words]}}
  KEYWORDS = {
    'eQTL': {
      'META': set(['expression', 'eqtl']),
      'GENE_SYMBOL': set(['gene', 'symbol', 'sym', 'genesym', 'genesymbol']),
      'ENTREZ_GENE_ID': set(['entrez', 'entrezid']),
      'ENSEMBL_ID': set(['ensembl', 'ensemblid']),
      'REFSEQ_ACC': set(['refseq', 'refseqacc', 'refseqaccession']),
      'GENBANK_ACC': set(['gb', 'acc', 'genbank', 'accession', 'genbankaccession'])},
    'SNP': {
      'META': set(['SNP', 'nucleotide']),
      'SNP_ID': set(['snp', 'id', 'rs', 'snpid', 'ncbi']),
      'CHROMOSOME': set(['chromosome', 'chrom', 'chr', 'ch']),
      'LOCATION': set(['mapinfo', 'map', 'info', 'loci', 'locus', 'loc', 'location', 'pos', 'position'])}
    }
  # Preferred order of gene row identification symbol
  EQTL_GENE_NAME_LIST = [
    'GENE_SYMBOL',
    'ENTREZ_GENE_ID',
    'ENSEMBL_ID',
    'REFSEQ_ACC',
    'GENBANK_ACC',
    ]
  COL_TYPE_RX = {
    # see: http://www.genenames.org/guidelines.html
    'GENE_SYMBOL': '[A-Z][a-zA-Z0-9-]+',
    'ENTREZ_GENE_ID': '\d+',
    # http://useast.ensembl.org/Help/View?id=143
    'ENSEMBL_ID': 'ENS[A-Z]{1,3}\d{11}',
    # http://www.ncbi.nlm.nih.gov/Sequin/acc.html    
    'GENBANK_ACC': '[A-Z]{1,5}[_-]?\d{2,8}',
    # http://www.ncbi.nlm.nih.gov/RefSeq/key.html
    'REFSEQ_ACC': '[A-Z]{2}_\d+',
    'SNP_ID': '(rs|cnvi)?\d+',
    'CHROMOSOME': '\d{1,2}',
    'LOCATION': '\d+',
  }
  for key in COL_TYPE_RX:
    COL_TYPE_RX[key] = re.compile(COL_TYPE_RX[key])

  def __init__(self, gpl_id, study_type=None, custom_parameters=None):
    """Initialize GPL.

    Args:
      gpl_id: str of GPL id like GPL\d+
      study_type: str in RECOGNIZED_STUDY_TYPES
      custom_parameters: {str: obj} of parameters specified by the user.
    """
    self.id = gpl_id
    self.populated = False
    self.loaded = False
    self.col_titles = []
    self.col_desc = {}
    self.row_desc = {}
    self.type = study_type
    self.special_cols = {}
    self.attrs = {}
    self.case_insensitive_row_id = {}
    # Set external parameters, update with custom user settings if they exist.
    self.parameters = GPL_SETTINGS.get(self.id, {}).copy()
    if custom_parameters:
      self.parameters.update(custom_parameters)

    # Populate self with meta data
    self._populate()

  def _check_id(self, line):
    """Raise error if platform id parsed from `line` does not match self.id."""
    try:
      platform_id = self.RX_PLATFORM.match(line).group(1)
    except:
      raise MalformedDataError, \
          "Cannot recognize GPL ID in line %s parsing %s" % (line, self.id)
    if platform_id != self.id:
      raise MalformedDataError, \
          "GPL ID %s differs from requested ID %s." % (platform_id, self.id)
  
  def _parse_brief(self, fp):
    """Parse GEO QUICK SOFT report for this GPL.

    Args:
      fp: [*str] of open line iterator to quick view, text form GPL description
    """
    # 1. Consume and check GPL ID
    # ==========
    line = fp.next().strip()
    self._check_id(line)

    # 2. Collect column title descriptions and attributes
    # ==========
    for line in fp:
      line = line.strip()

      # Is this a column description?
      m = self.RX_HEADER.match(line)
      if m:
        key, value = m.groups()
        self.col_desc[key] = value
        continue

      # Is this an attribute?
      m = self.RX_ATTR.match(line)
      if m:
        key, value = m.groups()
        self.attrs.setdefault(key, []).append(value)
        continue

      # Is this the end line? Exit loop.
      if line == self.HEAD_END_LINE:
        break

      # This line should never be reached in this loop.
      raise MalformedDataError, \
        "Unrecognized line '%s' while parsing %s" % (line, self)      

    # 3. Collect column titles from first row of data.
    # ==========
    line = fp.next()
    self.col_titles = line.strip().split("\t")

    # 4. Verify that all column titles have a description
    # ==========
    if len(self.col_titles) != len(self.col_desc):
      Log.warning("%d col_titles != %d col_desc for %s" % \
        (len(self.col_titles), len(self.col_desc), self))

  def __repr__(self):
    #TODO: include type
    return "[GPL %s (%d)]" % (self.id, id(self))

  def _populate(self):
    """Populate self with meta data, but not row descriptions."""

    # 1. Open GEO and load description of this GPL
    # ==========
    url = self.PTN_GPL_QUICK % {'id': self.id}
    handle = Download(url)
    http_fp = handle.read()
    self._parse_brief(http_fp)
    http_fp.close()

    # 2. Determine GPL type.
    # ==========
    # First try to predict what study type this is from meta data
    guessed_type = self._guess_type()
    # if given a type, check that this type matches guess.
    if self.type is not None:
      if guessed_type != self.type:
        Log.warning("Guessed type '%s' != given type '%s' for %s" % \
          (guessed_type, self.type, self))
      else:
        Log.info("Guessed type '%s' matches given type '%s' for %s" % \
          (guessed_type, self.type, self))
    else:
      self.type = guessed_type
      Log.info("Set type to guessed type '%s' for %s" % \
        (guessed_type, self))
    
    # 3. Identify special column titles (like gene symbol) for this GPL type
    # ==========
    special_cols = self.KEYWORDS[self.type]
    for name in special_cols:

      # Ignore "META" (it is only used for matching meta data, not col titles.
      if name == "META":
        continue
      
      # Override for custom-set paramaters.
      if name in self.parameters:
        if self.parameters[name] not in self.col_titles:
          raise MalformedDataError, \
            "Custom parameter %s=%s not for %s, col_titles=%s" % \
            (name, self.parameters[name], self, self.col_titles)
        mapped_name = self.parameters[name]
        
      # Automatically detect the best matching column title.
      else:
        mapped_name = self._find_column(name)

      self.special_cols[name] = mapped_name

    # Flag self as fully populated.
    self.populated = True

  def _guess_type(self):
    """Return a guess of the study type of this GPL from its description.

    Returns:
      str in RECOGNIZED_STUDY_TYPES
    """
    Log.info("Attempting to guess study type of GPL %s." % self)
    # Warn if no column title
    if not len(self.col_desc) > 1:
      raise NotPopulatedError, "Parse brief for %s before guessing type" % self

    # Match best type using keyword matching.
    type_ranks = {}
    for study_type in self.KEYWORDS:
      type_score = 0
      # Find keywords in column titles and their descriptions
      for col_title, keywords in self.KEYWORDS[study_type].items():
        for title, desc in self.col_desc.items():
          type_score += keyword_score(title, desc, keywords)
      # Also check the GPL attributes for meta keywords.
      keywords = self.KEYWORDS[study_type]['META']
      for key, values in self.attrs.items():
        type_score += 6*keyword_score(key, " ".join(values), keywords)
      # Update rankings
      type_ranks[study_type] = type_score

    # Guess top ranked study type. Log scores. Return best score.
    top = sorted(type_ranks, key=lambda x: type_ranks[x], reverse=True)[0]
    Log.info("Guessed type=%s for %s. Rankings: %s" % \
             (top, self, type_ranks))
    return top
    
  def _find_column(self, name):
    """Return column title most likely to represent the named column.
    GPL description headers must have been parsed to call this function.

    Args:
      name: str of column name in KEYWORDS for this GPL type.
    Returns:
      str of GPL column title name best matching `name` or None if no match
    """
    # Get static list of keywords.
    keywords = self.KEYWORDS[self.type][name]
    
    # Compute keyword match scores for all columns
    best_score = (None, 0)
    for title, desc in self.col_desc.items():
      score = keyword_score(title, desc, keywords)
      if score > best_score[1]:
        best_score = (title, score)
        
    # Warn if no good column has been found, else report mapping.
    title, score = best_score
    if title:
      desc = self.col_desc[title]
    else:
      desc = None
    if best_score[1] <= 1:
      # Do not return weak results.
      Log.warning("%s %s column not found with confidence. Best: '%s: %s'" % \
                  (self, name, title, desc))
      title = None
    else:
      Log.info("Mapped '%s' to col title '%s: %s' for %s. Confidence: %d" % \
               (name, title, desc, self, score))
                  
    # Return name of highest scoring GPL column title.
    return title

  def get_column(self, row_id, name):
    """Return value of special column at row id.

    Args:
      row_id: str of row ID in self.row_desc
      name: str of special row name in KEYWORDS[self.type]
    Returns:
      str of value in column mapped by `name` or None
    """
    if not self.populated:
      raise NotPopulatedError, "%s Row definitions not yet populated." % self

    # Verify that this special column exists for this GPL type.
    try:
      key = self.special_cols[name]
    except KeyError:
      raise KeyError, \
        "Special column %s not defined for %s of type %s. Defined cols: %s" % \
        (name, self, self.type, self.KEYWORDS[self.type].keys())

    try:
      row = self.row_desc[row_id]
    except KeyError:
      row_id = self.case_insensitive_row_id[row_id.lower()]
      row = self.row_desc[row_id]
    # Attempt to return this column value at this row, else return None
    try:
      value = row[key]
    except KeyError:
      value = None
    return value

  #TODO: populate from 
  
  def load(self):
    """Fetch GPL row definitions, load values into this object.
    """
    Log.info("Loading %s" % self)
    # Do not reload a populated GEO object.
    if self.loaded:
      Log.warning("%s already loaded." % self)
      return
    
    url = self.PTN_GPL % {'id': self.id}
    handle = Download(url)
    http_fp = handle.read()
    self._parse(http_fp)
    Log.info("Fetched %s while loading %s." % (url, self))
    http_fp.close()
    # Verify that at least one row description has been loaded.
    if len(self.row_desc) < 1:
      Log.warning("No row descriptions loaded for %s." % self)
    else:
      Log.info("Loaded %d row descriptions for %s." % (len(self.row_desc), self))
    self.loaded = True

  def _parse(self, fp):
    """Parse a GPL text representation.

    Args:
      fp: iter=>str of "#header" format lines
    """
    # 0. Consume and check GPL ID
    line = fp.next().strip()
    try:
      platform_id = self.RX_PLATFORM.match(line).group(1)
    except:
      raise MalformedDataError, \
          "Cannot recognize GPL ID in line %s parsing %s" % (line, self.id)
    if platform_id != self.id:
      raise MalformedDataError, \
          "GPL ID %s differs from requested ID %s." % (platform_id, self.id)
    
    # 1. Collect attributes and column titles.
    for line in fp:
      line = line.strip()
      m = self.RX_HEADER.match(line)
      # Only handle header lines in this loop.
      if not m:
        # This must be the end header line. 
        if line != self.HEAD_END_LINE:
          raise MalformedDataError, \
            "Unrecognized line '%s' while parsing %s" % (line, self.id)
        break
      key, value = m.groups()
      self.col_desc[key] = value

    # 2. Collect column titles from first row of data.
    line = fp.next()
    self.col_titles = line.strip().split("\t")
      
    # 3. 
    for line in fp:
      
      # Only strip end-of-line characters to avoid column misalignment.
      line = line.rstrip('\n\r')
      if line == self.TABLE_END_LINE:
        break
      
      # Split line by tabs.
      row = line.split('\t')
      
      # Create new dictionary for this row_id.
      row_id = row[0]
      self.row_desc[row_id] = {}

      # Populate the attributes for this row.
      for i in range(1, len(row)):
        value = row[i].strip()
        name = self.col_titles[i]
        # Ignore empty values.
        if value != "":
          self.row_desc[row_id][name] = value
          
      # Add row_id to case-insensitive map.
      if row_id.lower() in self.case_insensitive_row_id:
        Log.warning("Multiple case-insensitive row_ids map row_id %s for %s" %\
          (row_id.lower(), self))
      else:
        self.case_insensitive_row_id[row_id.lower()] = row_id
          
    # 4. Clean Up
    Log.info("Populated %s with %d row descriptions of %d columns." % \
             (self, len(self.row_desc), len(self.col_titles)))
    fp.close()

    
class GSM(object):
  """A GEO Sample definition.

  Attributes:
    id: GSM id like GSM23409
    populated: bool if self has values

    rx_title: regex with two groups: title and repetition
    subject: str of sample subject name
    rep: str of sample repetition identifier, if it exists
    attr: {str: [str]} of all associated "!Sample_*" GSM attribute values
  """
  
  def __init__(self, gsm_id, rx_title_str=None):
    if rx_title_str is None:
      rx_title_str = DFT_RX_TITLE_STR
    self.id = gsm_id
    self.populated = False
    self.rx_title = re.compile(rx_title_str)
    self.subject = None
    self.rep = None
    self.attr = {}

  def __repr__(self):
    return "[GSM %s (%d)]" % (self.id, id(self))

  def add_pair(self, key, value):
    """Populate self with key=>value as parsed from a Series Matrix file.

    Args:
      key: str
      value: str
    """
    self.attr.setdefault(key, []).append(value)
    if key == "title":
      m = self.rx_title.match(value)
      if m is None:
        Log.warning("Title '%s' of %s did not match regex %s." % \
                    (value, self, self.rx_title))
      else:
        self.subject, self.rep = m.groups()

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
      
