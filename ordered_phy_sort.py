#!/usr/bin/env python
import sys
import os
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("-input_path",
                    help="input directory path",
                    required=True)
parser.add_argument("-output_path",
                    help="output directory path" ,
                    required=True)

args = parser.parse_args()

input_path = args.input_path
output_path = args.output_path

if not os.path.isdir(output_path):
	os.makedirs(output_path)

phy_files = ("{}/{}".format(input_path,i) for i in os.listdir(input_path) if i.endswith(".phy"))
BASE_SET = set([" ", "\n", "A", "T", "G", "C", "?"])

order_list = ["and34", "m155", "fan66", "dasi523", "ranol61", "mas66",
"bema8", "mar55", "tafo69", "mery3", "jar3.7", "fia5.28",
"rano217", "pd37711", "pd37712", "pjx022", "pjx023",
"tandra4.17", "nofy610", "pt37709", "pt37710", "hih19",
"pc37498", "pc37516", "pc37706", "pc37707", "pc37708",
"pc37729", "bema23", "jam47", "kmetea73", "pv37713", "zah272",
"mm37780", "west", "4476", "6139", "7128", "oto30611"]

order_set = set(order_list)

keyset = set()
for phy_file in phy_files:
	lines = []
	sortkeys = {}

	with open(phy_file) as input_handle:
		for num, line in enumerate(input_handle):
			if num == 0:
				header = line
				header = [i for i in line.split(" ") if i]
			else:
				reformatted_line = line.split(" ", 1)
				if len(reformatted_line) > 1:
					new_seq=""
					for char in reformatted_line[1].upper():
						if char not in BASE_SET:
							char = "N"
						new_seq += char
					reformatted_line[1] = new_seq
				lines.append(reformatted_line)
				index = len(lines) - 1
				key = lines[index][0]
				if key in order_set:
					sortkeys[key] = index
				else:
					key = key.rsplit("_")[0]
					if key in order_set:
						lines[index][0] = key
						sortkeys[key] = index
				keyset.add(key)


	new_phy = output_path + "/" + phy_file.rsplit("/")[-1].rsplit(".phy")[0]+".phy"


	with open(new_phy, "w+") as out_handle:
		out_handle.write("   {}   {}".format(len(sortkeys), header[-1]))

		for item in order_list:
			if item in sortkeys:
				out_handle.write("{} {}".format(*lines[sortkeys[item]]))
