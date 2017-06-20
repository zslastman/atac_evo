#!/usr/bin/env python2.7
"""
Will tranform a given file name into a new file that incorporates the file path.
The output file will be produced in the same directory in which the program is run

usage: chage_file_name_to_path.py filePath
"""

import os
import sys
import shutil

inFile = sys.argv[1]
path_atoms = [k for k in inFile.split('/') if k != '.' and k != '..']
newFile = '_'.join(path_atoms)
shutil.copy(inFile, newFile)
