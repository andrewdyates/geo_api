#!/usr/bin/python
#
# Copyright 2012 Andrew D. Yates
#
"""Library and platform independent HTTP/FTP handles.

Cache directory is set locally by CACHE_DIR.
"""
import sys
import urllib
import urllib2
import ftplib
import re
import os

# patched, local version of gzip from Python 3 to handle http streams
import gzip3 as gzip

import download

# Hack to import application level Log object without adding it to global path
sys.path.append("..")
from logger import Log

# XXX: SET THIS TO YOUR LOCAL PATH
CACHE_DIR = "/Users/qq/biostat_data/cache"

def get_cache_name(url):
  """Return cache file name from url."""
  a = re.sub('^https?://','', url)
  a = re.sub('^ftp://', '', a)
  a = re.sub('[^a-zA-Z0-9.?]', "_", a)
  return a + ".cache"


class DownloadIter(object):
  """Iterator which reports download progress and caches downloads to file.

  Unlike normal file pointers, DownloadIter closes itself when its underlying
  HTTP buffer has been exhausted. This is done to ensure proper handling of 
  both the cache and the HTTP network connection.
  
  Subsequent calls to DownloadIter.close() after DownloadIter.closed = True have
  no effect. A premature call to DownloadIter.close() before it has closed itself
  will terminate the underlying HTTP connection and discard the open cache file.
  """
  # number of bytes to download before calling 'update' hook (100k)
  REPORT_SIZE = 131072
  def __init__(self, fp, size=None, cache=None, report=True, finalize=True, ftp=False):
    """Initialize self.

    Args:
      fp: open file pointer like object like HTTP connection
      size: int of expected bytes in download or None if unknown
      cache: str of cache file name or None if cache is disabled
      report: bool to report download status
      finalize: bool if close() is called before self.completed, finish 
        downloading before closing buffer and any cache
    """
    self.buffer = fp
    self.cache = cache
    self.size = size
    if self.size:
      # Convert to float to avoid integer division complications.
      self.size = float(self.size)
    self.report = report
    self.finalize = finalize

    self.bytes_read = 0
    self.bytes_reported = 0
    self.fp_out = None
    self.completed = False
    self.tmp_filepath = None
    self.dest_filepath = None
    self.closed = False

    if cache:
      # Finalized (completely downloaded) cache filename
      self.dest_filepath = os.path.join(CACHE_DIR, cache)
      # Add ".tmp" to end of cache filename to indicate cache is incomplete.
      self.tmp_filepath = self.dest_filepath + ".tmp"
      # Cache files are compressed (even if the underlying data is compressed)
      self.fp_out = gzip.open(self.tmp_filepath, "wb")

  def __repr__(self):
    return "[DownloadIter: cache=%s, %d bytes read, closed=%s, buffer=%s (%d)]" % \
      (self.cache, self.bytes_read, str(self.closed), self.buffer, id(self))

  def __iter__(self):
    """Call at start of iter read loops"""
    if self.completed or self.closed:
      Log.warning("Iterator opened on closed or completed %s" % self)
    return self

  def read(self, *args, **kwds):
    """Return string. Wrapper for direct access to underlying http buffer."""
    # The underlying buffer has been released. Just continue to return empty string
    #   as if the buffer were at EOF even though this handle has been "closed"
    if self.buffer is None:
      return ""
    block = self.buffer.read(*args, **kwds)
    self._handle_block(block)
    return block

  def _handle_block(self, block):
    """Internal method for caching and reporting download status.

    Args:
      block: str of data
    Returns:
      None if download is complete, else return unmodified, nonempty `block`
    """
    self.bytes_read += len(block)

    # Download complete: file.read() == "" means EOF. Return None.
    if block == "":
      self.completed = True
      if self.report:
        self._report()
      self.close()
      return None
    
    # Write block to cache (if enabled).
    if self.cache:
      self.fp_out.write(block)
    
    # Report download status (if enabled).
    if self.report:
      self.bytes_reported += len(block)
      if self.bytes_reported >= self.REPORT_SIZE:
        self.bytes_reported = 0
        self._report()
      
    return block
    
  def next(self):
    """Wrapper for next line iterator to underlying http buffer. 

    Returns:
      str of next line in http buffer
    """
    block = self.buffer.readline()
    if not self._handle_block(block):
      raise StopIteration
    return block

  def close(self):
    """Close any open file pointers, close and finalize cache file.
    """
    # Ignore repeated calls to close()
    if self.closed:
      Log.info("Redundant call to close(), Ignored for %s." % self)
      return
    else:
      Log.info("Closing %s..." % self)

    # Handle finalize requests to complete download to buffer.
    if self.finalize:
      if not self.completed and self.cache:
        Log.info("Finalizing download of %s." % self)
        # Read remaining buffer unconditionally. Use iterator if reporting.
        if self.report:
          while True:
            try:
              self.next()
            except StopIteration:
              break
        else:
          self.read()
        # If not closed in previous read(), try another read().
        if not self.closed:
          # This closes self since the previous read flushed the buffer.
          self.read()
        if not self.closed:
          Log.warning("Close sequence not completed as expected for %s." % self)
        # Exit: prior reads in the finalize process already closed self.
        return

    # self.buffer.close() causes bugs with FTP. Python sockets clean up after 
    #   themselves in garbage collection, so to remove the reference to buffer
    # self.buffer.close()
    self.buffer = None
    self.fp_out.close()

    if self.completed:
      Log.info("Download complete. %d bytes read." % (self.bytes_read))
      # Finalize cache.
      if self.cache:
        os.rename(self.tmp_filepath, self.dest_filepath)
        Log.info("Cache finalized as '%s'." % (self.dest_filepath))
    else:
      Log.info("Download closed before completion. %d bytes read." % \
               (self.bytes_read))
      # Flush cache.
      if self.cache:
        os.remove(self.tmp_filepath)
        Log.info("Incomplete cache '%s' deleted." % (self.tmp_filepath))
        
    # Flag self as closed to prevent redundant .close() calls.
    self.closed = True

  def _report(self):
    """Hook for reporting download status.

    Args:
      complete: bool download completed
    """
    if self.completed:
      Log.console("Download complete. %d bytes read." % (self.bytes_read))
    else:
      if self.size:
        percent_complete = self.bytes_read / self.size * 100
        if percent_complete >= 100:
          percent_complete = 99.99
        Log.console("%.2f%% downloaded. %d of %d bytes." % \
                    (percent_complete, self.bytes_read, self.size))
      else:
        Log.console("downloading: %d bytes so far..." % \
                    (self.bytes_read))


class CachedDownload(download.Download):
  """HTTP Download handler with automatic disk cache and compression support.
  Constructs and returns new DownloadIter fp from self.read().
  
  TODO: extend to support custom HTTP headers, other HTTP methods, etc.

  Attributes:
    status: int of response status
    headers: dict of response headers
    bytes: number of bytes downloaded
    expected_size: int of expected number of bytes to be downloaded
    url_fetched: str of url actually fetched (e.g., after following redirects)
  """
  # default headers for HTTP
  HEADERS = {
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
    'Accept-Charset': 'ISO-8859-1,utf-8;q=0.7,*;q=0.3',
    'Accept-Encoding': 'gzip',
    'Accept-Language': 'en-US,en;q=0.8'
  }
  RX_FTP = re.compile("^ftp://", re.I)
  
  def __init__(self, url, 
               req_data=None, req_headers=None, expected_size=None, report_status=True,
               write_cache=True, read_cache=True, finalize=True):
    """Specify HTTP request; init.

    Args:
      url: str of url to fetch
      req_data: {str:obj} of data to upload
      req_headers: {str:str} of additional HTTP headers
      expected_size: int of expected file size in bytes
      report_status: bool if to print download status to console
      write_cache: bool to write download to cache (if from network)
      read_cache: bool to read from cache (if it exists)
      finalize: bool if to finalize Download if it's closed before completion
    """
    # ftp or http?
    if self.RX_FTP.match(url):
      self.type = "ftp"
    else:
      self.type = "http"
    self.url = url
    self.bytes = 0
    self.req_data = req_data
    self.req_headers = req_headers
    self.report_status = report_status
    self.write_cache = write_cache
    self.read_cache = read_cache
    self.finalize = finalize
    
    if expected_size:
      self.expected_size = expected_size
    else:
      self.expected_size = None
      
    self.rsp_url = None
    self.headers = None
    self.status = None

    self.cache_name = get_cache_name(self.url)

  def __repr__(self):
    return "[CachedDownload %s (%d)]" % (self.url, id(self))

  def _fetch_from_cache(self):
    """Attempt to fetch a file from cache or return None.

    Returns:
      *iter of open, decompressed file pointer from cache or None.
    """
    filepath = os.path.join(CACHE_DIR, self.cache_name)
    if os.path.exists(filepath):
      return gzip.open(filepath, "rb")
    else:
      return None

  def fetch(self, data=None, headers=None):
    """Fetch http file from network.

    Args:
      headers: {str:str} of additional request HTTP headers
      data: {str:*} of data to be sent via HTTP
    Returns:
      [*str] of file pointer-like HTTP stream.
    """
    # Fetch request.
    if self.type == "http":
      rsp = self._fetch_http(data, headers)
    elif self.type == "ftp":
      rsp = self._fetch_ftp()
    else:
      Log.warning("Unknown type, cannot fetch %s for %s." % self.url, self)
      return None

    self.status = 200
    # Convert header keys into all lower case.
    self.headers = {}
    for key, value in dict(rsp.info()).items():
      self.headers[key.lower()] = value
    self.url_rsp = rsp.geturl()
    
    return rsp

  def _fetch_ftp(self):
    """Fetch from FTP. This may use the FTP library directly someday."""
    return urllib2.urlopen(self.url)

  def _fetch_http(self, data=None, headers=None):
    """Fetch from HTTP."""
    head = self.HEADERS
    if headers:
      head.update(headers)
    req = urllib2.Request(self.url, headers=head)
    if data:
      req.add_data(urllib.urlencode(data))

    try:
      rsp = urllib2.urlopen(req)
    except urllib2.URLError:
      self.status = rsp.code
      raise
    return rsp

  def read(self):
    """Return a file-pointer-like object to this resource.
    
    Returns:
      iter: file-pointer-like str line iterator (uncompressed)
    """
    # Attempt to retreive from cache if possible.
    if self.read_cache:
      fp = self._fetch_from_cache()
    else:
      fp = None
    if fp:
      Log.info("Fetched %s from cache." % self.url)
      return fp
    else:
      Log.info("Downloading %s from network." % self.url)
    
    # From HTTP, Fetch request and populate self with response.
    http_fp = self.fetch()
    # If compressed, wrap http handle in a gzip decompressor.
    if self.headers and "content-encoding" in self.headers and \
        self.headers["content-encoding"] == "gzip":
      zip_fp = gzip.GzipFile(fileobj=http_fp)
      fp = zip_fp
    else:
      fp = http_fp
      
    # Return download iterator from decompressed HTTP handle.
    if self.write_cache:
      cache = self.cache_name
    else:
      cache = None

    # Get expected download size in bytes.
    if self.headers and 'content-length' in self.headers:
      try:
        size = int(self.headers['content-length'])
      except:
        size = None
    else:
      size = None
      
    return DownloadIter( \
      fp, cache=cache, size=size, report=self.report_status, finalize=self.finalize)
