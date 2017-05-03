#!/usr/bin/env python3
from __future__ import print_function

import pprint
import sys
import urllib.parse
import urllib.request

from sphinx.ext.intersphinx import read_inventory_v1, read_inventory_v2


PACKAGES = {'python': 'http://docs.python.org/2',
            'numpy': 'http://docs.scipy.org/doc/numpy/',
            'scipy': 'http://docs.scipy.org/doc/scipy/reference/'}


def get_inventory(url):
    inv_url = urllib.parse.urljoin(url, 'objects.inv')
    with urllib.request.urlopen(inv_url) as f:
        line = f.readline().rstrip().decode('utf-8')
        if line == '# Sphinx inventory version 1':
            invdata = read_inventory_v1(f, url, urllib.parse.urljoin)
        elif line == '# Sphinx inventory version 2':
            invdata = read_inventory_v2(f, url, urllib.parse.urljoin)
        else:
            raise ValueError(line)
        return invdata


def main():
    if len(sys.argv) < 2:
        print('usage:', sys.argv[0], 'PACKAGE-NAME', file=sys.stderr)
        raise SystemExit(1)
    package = sys.argv[1]
    try:
        url = PACKAGES[package]
    except KeyError:
        print('unsupported package:', package, file=sys.stderr)
        raise SystemExit(1)
    result = get_inventory(url)
    pprint.pprint(result)


if __name__ == '__main__':
    main()
