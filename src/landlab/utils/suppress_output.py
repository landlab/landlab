import contextlib
import os


@contextlib.contextmanager
def suppress_output(out=True, err=True):
    """Suppress output from both stdout and stderr.

    Parameters
    ----------
    out : bool, optional
        Suppress stdout.
    err : bool, optional
        Suppress stderr.
    """
    null_fds, save_fds = {}, {}
    if out:
        null_fds["out"] = os.open(os.devnull, os.O_RDWR)
        save_fds["out"] = os.dup(1)
        os.dup2(null_fds["out"], 1)
    if err:
        null_fds["err"] = os.open(os.devnull, os.O_RDWR)
        save_fds["err"] = os.dup(2)
        os.dup2(null_fds["err"], 2)

    yield

    # Re-assign the real stdout/stderr back to (1) and (2)
    out and os.dup2(save_fds["out"], 1)
    err and os.dup2(save_fds["err"], 2)

    # Close the null files
    for fd in list(null_fds.values()) + list(save_fds.values()):
        os.close(fd)
