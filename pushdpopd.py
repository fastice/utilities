# grabbed this from https://trac.cnx.org/browser/trunk/python/pushd.py?rev=9428
#
import os

__dirstack = []
	
def pushd(directory=None):
    if directory:
        try:
            curdir = os.getcwd()
#            print(curdir)
        except OSError:
            # Weird, no current dir, let's hardcode a default
            curdir = "/tmp"
        __dirstack.append(curdir)
        os.chdir(directory)
    else:
        try:
            top = __dirstack.pop()
        except IndexError:
            print("pushd: No other directory.")
        else:
            __dirstack.append(os.getcwd())
            os.chdir(top)

	
def popd():
    try:
        os.chdir(__dirstack.pop())
    except IndexError:
        print("popd: Directory stack empty.")
            
