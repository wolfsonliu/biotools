#! /bin/env python3

# ------------------------------------------------------------------------------
# Package
import os, re
import pandas as pd
import urllib.request as url
import urllib.parse as urlparse

# ------------------------------------------------------------------------------
# Class

# ------------------
class KEGGAPI():
    '''
    KEGGAPI is used to abtain information from KEGG website.'''
    def __init__(self, operation, argument):
        self._url = urlparse.urljoin('http://rest.kegg.jp', operation, argument)
