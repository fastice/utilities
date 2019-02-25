# dols.py



from subprocess import check_output,STDOUT,CalledProcessError

def dols(lscmd):
    """ Excecute a unix ls command. You must specifcy ls in command. 
usage: dols("ls -d XYZ*")"""
    
    try:
        listing = check_output([lscmd],stderr=STDOUT,shell=True,executable='/bin/csh')
        returnCode = 0
    except CalledProcessError as ex:
        o=ex.output
        returnCode = ex.returncode
    if returnCode == 0 :
        listing=str(listing,'utf-8').split()
    else:
        listing=''

    return listing




