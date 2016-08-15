#!/usr/bin/env python

"""
.. module: redirector
   :synopsis: Redirect stdout or stderr

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import contextlib
import os
import sys

def fileno(file_or_fd):
    fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd

@contextlib.contextmanager
def stdout_redirected(to=os.devnull, stdout=None):
    """
    http://stackoverflow.com/a/22434262/190597 (J.F. Sebastian)
    """
    if to is None:
        yield stdout
        return
    if stdout is None:
       stdout = sys.stdout

    stdout_fd = fileno(stdout)
    # copy stdout_fd before it is overwritten
    #NOTE: `copied` is inheritable on Windows when duplicating a standard stream
    with os.fdopen(os.dup(stdout_fd), 'wb') as copied: 
        stdout.flush()  # flush library buffers that dup2 knows nothing about
        try:
            os.dup2(fileno(to), stdout_fd)  # $ exec >&to
        except ValueError:  # filename
            with open(to, 'wb') as to_file:
	              os.dup2(to_file.fileno(), stdout_fd)  # $ exec > to
        try:
            yield stdout # allow code to be run with the redirected stdout
        finally:
            # restore stdout to its previous value
            #NOTE: dup2 makes stdout_fd inheritable unconditionally
            stdout.flush()
            if stdout is not None:
                os.dup2(copied.fileno(), stdout_fd)  # $ exec >&copied

@contextlib.contextmanager
def stderr_redirected(to=os.devnull, stderr=None):
    """
    http://stackoverflow.com/a/22434262/190597 (J.F. Sebastian)
    """
    if to is None:
        yield stderr
        return

    if stderr is None:
       stderr = sys.stderr

    stderr_fd = fileno(stderr)
    # copy stderr_fd before it is overwritten
    #NOTE: `copied` is inheritable on Windows when duplicating a standard stream
    with os.fdopen(os.dup(stderr_fd), 'wb') as copied: 
        stderr.flush()  # flush library buffers that dup2 knows nothing about
        try:
            os.dup2(fileno(to), stderr_fd)  # $ exec >&to
        except ValueError:  # filename
            with open(to, 'wb') as to_file:
	              os.dup2(to_file.fileno(), stderr_fd)  # $ exec > to
        try:
            yield stderr # allow code to be run with the redirected stderr
        finally:
            # restore stderr to its previous value
            #NOTE: dup2 makes stderr_fd inheritable unconditionally
            stderr.flush()
            os.dup2(copied.fileno(), stderr_fd)  # $ exec >&copied

if __name__ == '__main__':
    print ( "test. this should appear." )
    with stdout_redirected(  ):
			print ( "redirected. this must not appear." )
    with stdout_redirected( to = None ):
			print ( "redirected trivially. this should appear." )
    print ( "test. this should appear." )
