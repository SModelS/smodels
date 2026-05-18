#!/usr/bin/env python3

from ptools.helpers import py_dump

with open("matrix.cov","rt") as f:
    t = eval(f.read())
    py_dump ( t, "new.cov", stop_at_level=1 )


