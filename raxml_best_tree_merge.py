#!/usr/bin/env python
import argparse
import os

"""This script was for merging raxml trees developed by Ryan Culligan and slightly modified"""

parser = argparse.ArgumentParser(description='')
parser.add_argument('-input_path', action="store", dest="input_path", required=True)
parser.add_argument('-output_path', action="store", dest="output_path",required=True)
args = parser.parse_args()
input_path = args.input_path
output_path = args.output_path

tree_string = ""
first = True
for directory in os.listdir(input_path):
    directory = os.path.join(input_path, directory)
    
    if os.path.isdir(directory):
        for best_tree in os.listdir(directory):
            if best_tree.startswith("RAxML_bestTree."):
                file_path = os.path.join(directory, best_tree)
                print(file_path)
                with open(file_path) as file_handle:
                    if first:
                        tree_string += file_handle.read()
                        print tree_string
                        first = False

                    else:
                        tree_string += "\n"+file_handle.read()


with open(output_path, "w+") as output_handle:
    output_handle.write(tree_string) 

print("WROTE FILE TO {}".format(args.output_path))
