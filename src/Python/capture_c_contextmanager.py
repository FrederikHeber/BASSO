# This is the answer to a question asked on stackoverflow:
# https://stackoverflow.com/questions/35745541/how-to-get-printed-output-from-ctypes-c-functions-into-jupyter-ipython-notebook
#
# This is the code by Denilson Sá Maia

import io
import os
import sys
import tempfile
from contextlib import contextmanager

from basso_capture import fflush_stderr, fflush_stdout, \
    disable_stderr_buffering, disable_stdout_buffering

@contextmanager
def capture_std_out_err(encoding='utf8'):
    # Flushing, it's a good practice.
    sys.stdout.flush()
    sys.stderr.flush()
    fflush_stdout()
    fflush_stderr()

    # We need to use a actual file because we need the file descriptor number.
    with tempfile.TemporaryFile(buffering=0) as temp:
        # Saving a copy of the original stdout.
        prev_sys_stdout = sys.stdout
        prev_stdout_fd = os.dup(1)
        os.close(1)

        prev_sys_stderr = sys.stderr
        prev_stderr_fd = os.dup(2)
        os.close(2)

        # Duplicating the temporary file fd into the stdout fd.
        # In other words, replacing the stdout.
        os.dup2(temp.fileno(), 1)
        os.dup2(temp.fileno(), 2)

        # Replacing sys.stdout for Python code.
        #
        # IPython Notebook version of sys.stdout is actually an
        # in-memory OutStream, so it does not have a file descriptor.
        # We need to replace sys.stdout so that interleaved Python
        # and C output gets captured in the correct order.
        #
        # We enable line_buffering to force a flush after each line.
        # And write_through to force all data to be passed through the
        # wrapper directly into the binary temporary file.
        temp_wrapper = io.TextIOWrapper(
            temp, encoding=encoding, line_buffering=True, write_through=True)
        sys.stdout = temp_wrapper
        sys.stderr = temp_wrapper

        # Disabling buffering of C stdout.
        disable_stdout_buffering()
        disable_stderr_buffering()

        yield

        # Must flush to clear the C library buffer.
        fflush_stdout()
        fflush_stderr()

        # Restoring stdout.
        os.dup2(prev_stdout_fd, 1)
        os.close(prev_stdout_fd)
        sys.stdout = prev_sys_stdout

        os.dup2(prev_stderr_fd, 2)
        os.close(prev_stderr_fd)
        sys.stderr = prev_sys_stderr

        # Printing the captured output.
        temp_wrapper.seek(0)
        print(temp_wrapper.read(), end='')
