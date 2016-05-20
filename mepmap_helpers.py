# Minet small helper functions

# Import modules
import sys

def s_out(string):
    sys.stdout.write(string)
    sys.stdout.flush()


def s_err(string):
    sys.stderr.write(string)
    sys.stderr.flush()


def chunks(lst,n):
    return [ lst[i::n] for i in range(n) ]

def joinit(iterable, delimiter):
    """Intersperse iterable with a delimiter"""
    it = iter(iterable)
    yield next(it)
    for x in it:
        yield delimiter
        yield x
