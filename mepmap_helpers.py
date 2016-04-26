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
