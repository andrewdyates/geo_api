#!/usr/bin/python
"""GZip handler which closes its underlying fileobj when closed. 
Uses local version of gzip3 library backported from Python 3.
"""
import gzip3 as gzip

class Gzipper(gzip.GzipFile):
  """Wrapper for gzip which closes underlying streams when closed."""
  def close(self, *args, **kwds):
    """Close underlying fileobj, then close self."""
    if self.fileobj:
      self.fileobj.close()
    return super(Gzipper, self).close(*args, **kwds)
    
