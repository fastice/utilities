import os
import time
import sys
from datetime import datetime
import utilities as u


def runMyThreads(threads, maxThreads, message, delay=0.2, prompt=False):
    """ Loop through a list of  -- threads -- starting each one.
    Allow -- maxThreads -- running at once.
    In each loop iteration check status and print number running along with
    -- message.
    """
    #
    # Optional prompt
    #
    if prompt:
        while(1):
            ans = input(
                f'\n\033[1mRun these {str(len(threads))} jobs [y/n] \033[0m\n')
            if ans.lower() == 'y':
                break
            if ans.lower() == 'n':
                u.myerror("User prompted abort")
    notDone = True
    nRun = count = 0
    running = []
    # make sure always calling from home directory
    home = os.getcwd()
    # delay
    if delay < 0.2:
        delay = 0.2
    # format codes
    bs = '\033[1m'
    norm = '\033[0m'
    grs = '\033[1;42m'
    bls = '\033[1;44m'
    # time for counter
    start = datetime.now()
    # loop counter
    n = 0
    while notDone:
        #
        # Start a thread is < maxThreads
        #
        if n % 1 == 0:
            timeElapsed = datetime.now() - start
            print(grs, message, '(', maxThreads, ')', norm, ': nRunning ', bs,
                  f'{nRun:5}', norm, ' nStarted ', bs, f'{count:5}',
                  norm, 'nToGo ', bs, f'{len(threads)-count:5}', '  ', bls,
                  timeElapsed, norm, '     '.ljust(maxThreads+8), end='\r')
            print(grs, message, '(', maxThreads, ')', norm, ': nRunning ', bs,
                  f'{nRun:5}', norm, 'nStarted ', bs, f'{count:5}',
                  norm, 'nToGo ', bs, f'{len(threads)-count:5}', '  ', bls,
                  timeElapsed, norm, '     ', end='')
            sys.stdout.flush()
        #
        # Run as many threads as can be started
        while count < len(threads) and nRun < maxThreads:
            # run thread and always make sure to return to current directory
            os.chdir(home)
            print('.', end='')
            sys.stdout.flush()
            threads[count].start()
            running.append(threads[count])
            nRun += 1
            count += 1
            time.sleep(delay)
            #
            # check status of threads
            #
        toRemove = []
        #
        # loop through running thread to find threads that are done
        for t in running:
            if not t.isAlive():
                toRemove.append(t)
        #
        # update list of running
        for t in toRemove:
            nRun -= 1
            running.remove(t)
        time.sleep(1)
        #
        print('', end='\r')
        n += 1
        if nRun == 0 and count >= len(threads):
            notDone = False
    print('--\n')
    u.myalert('Threads Done')
    return
