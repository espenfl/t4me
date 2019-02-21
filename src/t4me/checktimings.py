#!/usr/bin/python
"""Contains some information needed by timing function."""
# pylint: disable=invalid-name

import pstats
p = pstats.Stats('timings')
# keys available for sort stats
# 'calls'	call count
# 'cumulative'	cumulative time
# 'cumtime'	cumulative time
# 'file'	file name
# 'filename'	file name
# 'module'	file name
# 'ncalls'	call count
# 'pcalls'	primitive call count
# 'line'	line number
# 'name'	function name
# 'nfl'	        name/file/line
# 'stdname'	standard name
# 'time'	internal time
# 'tottime'	internal time

# only print the first 20 lines
p.strip_dirs().sort_stats('cumtime').print_stats(20)
