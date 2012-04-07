#!/usr/bin/python
#
# Copyright 2012 Andrew D. Yates
#
"""Library and platform independent HTTP/FTP handles.

LocalDownload: download to PC using ordinary HTTP calls and disk writes.

Notes:

The gzip file handle cannot directly accept an http "file handle" without raising an exception like:

  File "/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/gzip.py", line 252, in read
    self._read(readsize)
  File "/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/gzip.py", line 279, in _read
    pos = self.fileobj.tell()   # Save current position
AttributeError: addinfourl instance has no attribute 'tell'

The fix is to the the backported gzip library from Python 3.
"""
import sys
import urllib
import urllib2
import re
import logging

# patched, local version of gzip from Python 3 to handle http streams
import gzip3 as gzip

import download


# class DownloadIter:
#   """Iterator which reports download progress."""
#   pass


class LocalDownload(download.Download):
  """Default download object for local machine use. 

  Attributes:
    status: int of response status
    headers: dict of response headers
    bytes: number of bytes downloaded
    expected_size: int of expected number of bytes to be downloaded
    url_fetched: str of url actually fetched (e.g., after following redirects)
  """
  # number of bytes to download before calling 'update' hook (100k)
  REPORT_BYTES = 131072
  # default headers for HTTP
  HEADERS = {
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
    'Accept-Charset': 'ISO-8859-1,utf-8;q=0.7,*;q=0.3',
    'Accept-Encoding': 'gzip',
    'Accept-Language': 'en-US,en;q=0.8'
  }
  RX_FTP = re.compile("^ftp://", re.I)
  
  def __init__(self, url, req_data=None, req_headers=None, expected_size=None):
    """Specify HTTP request; init.

    Args:
      url: str of url to fetch
      req_data: {str:obj} of data to upload
      req_headers: {str:str} of additional HTTP headers
      expected_size: int of expected file size in bytes
      **kwds: consume extra keyword arguments (for compatibility)
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
    if expected_size:
      self.expected_size = expected_size
    else:
      self.expected_size = None
      
    self.rsp_url = None
    self.headers = None
    self.status = None
    

  def fetch(self, data=None, headers=None):
    """Return file pointer like object of fetched data.

    Args:
      headers: {str:str} of additional request HTTP headers
      data: {str:*} of data to be sent via HTTP
      **kwds: consume extra keyword arguments (for compatibility)
    """
    # don't set default headers for FTP; I don't know 
    if self.type != "ftp":
      head = self.HEADERS
    else:
      head = {}
    if headers:
      head.update(headers)
    if head:
      req = urllib2.Request(self.url, headers=head)
    else:
      req = urllib2.Request(self.url)
    if data:
      req.add_data(urllib.urlencode(data))

    # Fetch request.
    try:
      rsp = urllib2.urlopen(self.url)
    except urllib2.URLError:
      self.status = rsp.code
      raise
    self.status = 200
    self.headers = dict(rsp.info())
    self.url_rsp = rsp.geturl()
    return rsp


  def read(self):
    """Return a file-pointer-like object to this http document.

    Args:
      **kwds: arguments passed to self.fetch
    Returns:
      iter: file-pointer-like str line iterator (uncompressed)
    """
    # Fetch request and populate self with response.
    http_fp = self.fetch()
    # If compressed, decompress string before return
    if "Content-Encoding" in self.headers and \
        self.headers["Content-Encoding"] == "gzip":
      zip_fp = gzip.GzipFile(fileobj=http_fp)
      return zip_fp
    else:
      return http_fp

  def _report(self, complete=False):
    """Hook for reporting download status.

    Args:
      complete: bool download completed
    """
    if complete:
      logging.info("%s download complete. %d bytes." % (self.url, self.bytes))
    else:
      if self.expected_size:
        percent_complete = self.bytes / float(self.expected_size) * 100
        if percent_complete >= 100:
          percent_complete = 99.99
          logging.info("%s %.2f%% downloaded: %d of %d bytes." % \
            (self.url, percent_complete, self.bytes, self.expected_size))
        else:
          logging.info("%s downloading: %d bytes so far..." % \
            (self.url, self.bytes))

  def save(self, dest):
    """Save HTTP data to a destination (e.g., file, blobstore) with reporting.
    Intended for large downloads.

    Args:
      dest: str of destination identifier
      **kwds: arguments passed to self.fetch
    """
    http_fp = self.fetch()
    out_fp = open(dest, "wb")

    self.bytes = 0
    n = 0
    for blob in http_fp:
      out_fp.write(blob)
      self.bytes += len(blob)
      n += len(blob)
      # call report hook if enough bytes have been downloaded
      if n >= self.REPORT_BYTES:
        self._report()
        n = 0
    # Download complete. Report status.
    self._report(complete=True)

    out_fp.close()
