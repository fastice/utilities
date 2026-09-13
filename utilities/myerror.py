# dols.py
import sys
def myerror(message,myLogger=None):
    """ print error and exit with a FAILING status (1).

    Was a bare sys.exit(), which exits 0 -- so every fatal error here looked
    like success to subprocess.run().returncode, csh $status and && chains,
    and silently vanished from cron summaries.

    Note: inside a worker thread this raises SystemExit in that thread only;
    the process exit code is unaffected either way.
    """
    print('\n\t\033[1;31m *** ',message,' *** \033[0m\n')
    if myLogger != None :
        myLogger.logError(message)
    sys.exit(1)




