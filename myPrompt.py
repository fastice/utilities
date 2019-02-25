# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 07:34:44 2018

@author: ian
"""
from utilities import myerror
def myPrompt(promptText,abort=None) :
    ''' Prompt for y or n and return true or false - optionally abort if abort=True'''
    while(1) :
        ans=input('\n\033[1m {0:s}  [y/n] \033[0m\n'.format(promptText))
        if ans.lower() == 'y'  :
            return True
        if ans.lower() == 'n'  :
            if abort :
                myerror("User prompted abort")
            return False