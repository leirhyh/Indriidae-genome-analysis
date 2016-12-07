#!/usr/bin/env/python
#!/opt/anaconda1anaconda2anaconda3/bin/python
# -*- coding: utf-8 -*-

"""
(c) 2015 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 07 August 2012 21:08 PDT (-0700)
"""


import os
import glob
import argparse
import multiprocessing
from Bio import AlignIO
from collections import Counter
from phyluce.helpers import is_dir, FullPaths, get_file_extensions

#import pdb#!/opt/anaconda1anaconda2anaconda3/bin/python
# -*- coding: utf-8 -*-

"""
(c) 2015 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 07 August 2012 21:08 PDT (-0700)
"""


import os
import glob
import argparse
import multiprocessing
from Bio import AlignIO
from collections import Counter
from phyluce.helpers import is_dir, FullPaths, get_file_extensions

#import pdb
import shutil

def get_args():
    parser = argparse.ArgumentParser(
            description="""Count the number of informative sites in a given set of alignment"""
        )
    parser.add_argument(
            '--input',
            type=str,
            dest="input",
            default=None,
            help="""The input directory name""",
            required=True
        )
    parser.add_argument(
            '--output',
            type=str,
            default=None,
            help="""The output filename""",
            required=True
        )
    parser.add_argument(
            "--input-format",
            dest="input_format",
            choices=['fasta', 'nexus', 'phylip', 'clustal', 'emboss', 'stockholm'],
            help="""The input alignment format""",
            default='phylip'
        )
    parser.add_argument(
            "--cores",
            type=int,
            default=1,
            help="""The number of cores to use.""",
        )
    parser.add_argument(
            "--min_all_seqs_parsimony",
            type=int,
            default=25,
            help="""The maximum percentage of missing parsimony sites for any individual, after which file isn't printed""")

    parser.add_argument(
            "--output_phy_dir",
            default="phy_filtered_by_parsimony_default",
            help="""passing phylip files go here""")

    parser.add_argument(
            "--informative_sites_path",
            default=None,
            type=str,
            help="""path""")

    parser.add_argument(
            "--output_stats_csv",
            default="missing_sites_percentage_default.csv",
            help="""passing phylip files go here""")
    parser.add_argument(
            "--output_id_sites",
            default="percent_missing_parsimony_sites_default.csv",
            help="""passing phylip files go here""")


    return parser.parse_args()


def get_files(input_dir, input_format):
    alignments = []
    for ftype in get_file_extensions(input_format):
        alignments.extend(glob.glob(os.path.join(input_dir, "*{}".format(ftype))))
    return alignments


def get_informative_sites(count):
    counter = 0
    for i in count:
        if i == '?':
            counter += 1
    percent_missing = (counter / float(len(count))) * 100
    if percent_missing == 0:
        percent_missing = True

    # remove gaps
    del count['-']
    # remove N
    del count['N']
    # remove ?
    del count['?']

    sufficient_sites = len(count)
    if sufficient_sites >= 2:
        sufficient_sequences = sum([1 for i in count.values() if i >= 2])
        if sufficient_sequences >= 2:
            return percent_missing
    return False

def get_differences(count):
    # remove gaps
    del count['-']
    # remove N
    del count['N']
    # remove ?
    del count['?']
    # remove X
    del count['X']
    sufficient_sites = len(count)
    # counted, different = (1,1)
    if sufficient_sites >= 2:
        return (1, 1)
    # counted, not different = (1,0)
    elif sufficient_sites >= 1 and count.most_common()[0][1] > 1:
        return (1, 0)
    # not counted, not different = (0,0)
    else:
        return (0, 0)

def worker(work):
    args, f = work
    aln = AlignIO.read(f, args.input_format)
    name = os.path.basename(f)
    informative_sites = []
    percent_question_mark_at_site = []
    differences = []
    counted_sites =  []

    alignment_dict = {}
    all_alignment_dict = {}

    for index in aln:
        alignment_dict[str(index.id)] = []
        all_alignment_dict[str(index.id)] = []

    for idx in xrange(aln.get_alignment_length()):
        col = aln[:, idx].upper()
        count = Counter(col)
        percent_missing = get_informative_sites(count)
        for index in aln:
            all_alignment_dict[str(index.id)].append(index[idx])
        #print "NOT MISSING", percent_missing
        if percent_missing:
            informative_sites.append(1)
            for index in aln:
                alignment_dict[str(index.id)].append(index[idx])
            if percent_missing:
                percent_missing = 0
            percent_question_mark_at_site.append(percent_missing)
        else:
            informative_sites.append(0)
        diff = get_differences(count)
        if diff == (1, 1):
            counted_sites.append(1)
            differences.append(1)
        elif diff == (1, 0):
            differences.append(0)
            counted_sites.append(1)
        else:
            differences.append(0)
            counted_sites.append(0)

    if len(percent_question_mark_at_site):
        mean_percent_question_at_informative_site = str(percent_question_mark_at_site)
    else:
        mean_percent_question_at_informative_site = "undefined"

    bad_set = ["-", "?", "N"]
    new_alignment_dict = {}
    for individual in all_alignment_dict:
        if individual not in alignment_dict:
            new_alignment_dict[individual] = 0.00
        length_int = len(alignment_dict[individual])
        float_len = float(length_int)
        if length_int == 0:
            new_alignment_dict[individual] = 100
        else:
            new_alignment_dict[individual] = (len([i for i in alignment_dict[individual] if i in bad_set]) /  float_len) * 100
    new_alignment_dict["ALIGNMENT_FILE"] = f
    new_alignment_dict["ALIGNMENT"] = aln
    alignment_dict["ALIGNMENT_FILE"] = f
    alignment_dict["ALIGNMENT"] = aln
    for individual in all_alignment_dict:
        float_len = float(len(all_alignment_dict[individual]))
        if float_len == 0:
            all_alignment_dict[individual] = [100]
        else:
            all_alignment_dict[individual] = [len([i for i in all_alignment_dict[individual] if i in bad_set]) /  float_len * 100]
    all_alignment_dict["ALIGNMENT_FILE"] = f
    all_alignment_dict["ALIGNMENT"] = aln

    return (name, aln.get_alignment_length(), sum(informative_sites), sum(differences), sum(counted_sites), mean_percent_question_at_informative_site, new_alignment_dict, all_alignment_dict, alignment_dict)


def main():
    args = get_args()

    path = args.output_phy_dir
    new_path = path
    count = 0
    while os.path.isdir(new_path):
        count += 1
        new_path = "{}_{}".format(path, count)
    os.makedirs(new_path)
    args.output_phy_dir = new_path

    path = args.informative_sites_path
    new_path = path
    count = 0
    while os.path.isdir(new_path):
        count += 1
        new_path = "{}_{}".format(path, count)
    os.makedirs(new_path)
    args.informative_sites_path = new_path

    work = [(args, f) for f in get_files(args.input, args.input_format)]
    if args.cores <= 1:
        results = map(worker, work)
    elif args.cores > 1:
        pool = multiprocessing.Pool(args.cores)
        results = pool.map(worker, work)
    if args.output:
        outf = open(args.output, 'w')
        outf.write("locus,length,informative_sites,differences,counted-bases,mean_percent_question_at_informative_site\n")
    else:
        print "locus\tlength\tinformative_sites\tdifferences\tcounted-bases,mean_percent_question_at_informative_site"
    total_sites = []
    total_differences = []
    all_counted_sites = []
    mean_percent_question_at_informative_sites = []
    alignment_dicts = []
    all_alignment_dict = []
    seq_alignment_dict = []
    master_alignments = {}
    for locus in results:
        total_sites.append(locus[2])
        total_differences.append(locus[3])
        all_counted_sites.append(locus[4])
        mean_percent_question_at_informative_sites.append(locus[5])
        alignment_dicts.append(locus[6])
        all_alignment_dict.append(locus[7])
        seq_alignment_dict.append(locus[8])
        if not args.output:
            print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(locus[0], locus[1], locus[2], locus[3], locus[4], locus[5])
        else:
            outf.write("{0},{1},{2},{3},{4},{5}\n".format(locus[0], locus[1], locus[2], locus[3], locus[4], locus[5]))

    counter = 0

    print args.output_id_sites, "OUTPUT_ID_SITES"

    with open(args.output_id_sites, "w+") as output_handle:
        keys = []
        for num1, alignment_dict in enumerate(alignment_dicts):
            all_low = True
            if num1 == 0:
                output_handle.write("ALIGNMENT_FILE")
                for key in alignment_dict:
                    keys.append(key)
                    if key not in ("ALIGNMENT", "ALIGNMENT_FILE"):
                        output_handle.write(",{}".format(key))
                output_handle.write("\n")

            for num2, key in enumerate(keys):
                if not num2:
                    output_handle.write(alignment_dict["ALIGNMENT_FILE"])
                if key in ("ALIGNMENT", "ALIGNMENT_FILE"):
                    continue

                if key in alignment_dict:
                    value = alignment_dict[key]
                else:
                    value = 0.00

                str_value = ",{}".format(value)
                output_handle.write(str_value)

                if value > args.min_all_seqs_parsimony:
                    all_low = False

            output_handle.write("\n")

            if all_low:
                counter += 1
                print(alignment_dict["ALIGNMENT_FILE"])
                shutil.copy(alignment_dict["ALIGNMENT_FILE"], args.output_phy_dir)

    with open(args.output_stats_csv, "w+") as output_handle:
        for num, alignment_dict in enumerate(all_alignment_dict):
            all_low = True
            if num == 0:
                output_handle.write("ALIGNMENT_FILE")
                for key in alignment_dict:
                    if key not in ("ALIGNMENT", "ALIGNMENT_FILE"):
                        output_handle.write(",{}".format(key))
                output_handle.write("\n")
            for num2, key in enumerate(alignment_dict):
                if num2 == 0:
                    output_handle.write(alignment_dict["ALIGNMENT_FILE"])
                if key in ("ALIGNMENT", "ALIGNMENT_FILE"):
                    continue
                #print key
                value = ",{}".format(alignment_dict[key])
                output_handle.write(value)
            output_handle.write("\n")
    print 'in the', counter, "sequences above, all individuals have at least ", args.min_all_seqs_parsimony, "% of parsimony sites intact"

    print "\n\n"
    for alignment_dict in alignment_dicts:
        for key in alignment_dict:
            if key in ("ALIGNMENT", "ALIGNMENT_FILE"):
                continue
            if key not in master_alignments:
                master_alignments[key] = []
            master_alignments[key].append(alignment_dict[key])

    sequence_alignments = {}
    for alignment_dict in seq_alignment_dict:
        sequence_alignments[alignment_dict["ALIGNMENT_FILE"]] = {}
        for key in alignment_dict:

            if key in ("ALIGNMENT", "ALIGNMENT_FILE"):
                continue

            if key not in sequence_alignments:
                sequence_alignments[alignment_dict["ALIGNMENT_FILE"]][key] = []
            sequence_alignments[alignment_dict["ALIGNMENT_FILE"]][key] += alignment_dict[key]

    for individual in sequence_alignments:
        gene = sequence_alignments[individual]
        key_base = individual.rsplit("/")[-1].rsplit(".")[0]
        new_key = args.informative_sites_path + "/" + key_base + "_only_informative.phy"
        with open(new_key, "w+") as output_handle:
            for num, key in enumerate(gene):
                name = key_base
                if num == 0:
                    length = len(gene[key])
                    seqs = len(gene)
                    output_handle.write("{}  {}\n".format(seqs,length))
                output_handle.write("{}   {}\n".format(key, "".join(gene[key])))

    print "\n\n"

    print "Total sites = {0}; Sites per locus = {1:.2f}; Total differences = {2}; Differences per locus = {3:.2f}; All sites checked for differences = {4}".format(
        sum(total_sites),
        sum(total_sites) / float(len(total_sites)),
        sum(total_differences),
        sum(total_differences) / float(len(total_differences)),
        sum(all_counted_sites)
    )
    print "\n\n"

if __name__ == '__main__':
    main()
