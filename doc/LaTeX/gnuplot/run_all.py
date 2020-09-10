#!/usr/bin/env python3.3
import os
from os import system

rootDir = os.getcwd()

for directory in os.listdir(rootDir):
    if os.path.isdir(directory):
        currentDir = os.path.join(rootDir, directory)
        for files in os.listdir(currentDir):
            if files.endswith('.gp'):
                os.chdir(currentDir)
                command ='gnuplot %s' % files
                system(command)				
                os.chdir(rootDir)
                
                

for directory in os.listdir(rootDir):
    if os.path.isdir(directory):
        currentDir = os.path.join(rootDir, directory)
        print(currentDir)
        for files in os.listdir(currentDir):
            if files.endswith('.eps'):
                print(files)
                os.chdir(currentDir)
                command ='epstopdf %s' % files         
                system(command)
                fileName ='%s' % files  		
                #os.remove(fileName)
                fileNameBase, fileExt = os.path.splitext(fileName)
                #os.rename('%s.eps.pdf' % fileNameBase, '%s.pdf' % fileNameBase)
                os.chdir(rootDir)
                                                          