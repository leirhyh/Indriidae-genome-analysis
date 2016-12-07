#!/usr/bin/env python
import argparse
import os
import shutil

"""this python script to pull all the bootstrap files to one
 folder for Astral species tree bootstrap, developed with Ryan Culligan"""

parser = argparse.ArgumentParser(description='')   # standard input and output format using argparse
parser.add_argument('-input_path', action="store", dest="input_path", required=True)
parser.add_argument('-output_path', action="store", dest="output_path",required=True)
args = parser.parse_args()
input_path = args.input_path
output_path = args.output_path

path = os.getcwd()
outputPath = os.path.join(path, output_path)

if not os.path.isdir(output_path): # check if the folder is exist, if not, make an new folder.
	os.makedirs(output_path)

for directory in os.listdir(input_path):
    directory = os.path.join(input_path, directory)
    if os.path.isdir(directory):
        for best_tree in os.listdir(directory):
            if best_tree.startswith("RAxML_bootstrap."):
                file_path = os.path.join(directory, best_tree)
                
                shutil.copy2(file_path, outputPath )  # copy the single file to a specific directory


