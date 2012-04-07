"""Library and platform independent HTTP/FTP and GZip handles.

Customize this file to return the appropriate handlers per platform.
"""
#import local_download
import cached_download
import gzipper

#Download = local_download.LocalDownload
Download = cached_download.CachedDownload
Gzipper = gzipper.Gzipper
