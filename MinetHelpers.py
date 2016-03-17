# Minet small helper functions

# Import modules
import sys

def sWrite(string):
    sys.stdout.write(string)
    sys.stdout.flush()


def sError(string):
    sys.stderr.write(string)
    sys.stderr.flush()


def Chunks(lst,n):
    return [ lst[i::n] for i in range(n) ]
