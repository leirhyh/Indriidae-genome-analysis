#!/usr/bin/env python

import os
import subprocess as sp
import argparse
import uuid
import sys
import time 
import threading
from collections import deque

"""This script was for parallel running raxml of the phylip files by taxon developed by Ryan Culligan"""


try:
    import Queue
except ImportError:
    import queue as Queue

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class myThread (threading.Thread):
    def __init__(self, threadID, name, q):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.q = q
    def run(self):
        print("Starting " + self.name)
        process_data(self.name, self.q)#,self.alpha)
        print("Exiting " + self.name)

def process_data(threadName, q):
    while not exitFlag:
        queueLock.acquire()
        if not workQueue.empty():
            data = q.get()
            queueLock.release()
            for item in data:
                sp.call('{}'.format(item), shell=True)
        else:
            queueLock.release()
        time.sleep(1)


def remove_ignored(path, output_path, visited=set()):
    
    #make directory
    if path not in visited:
        print(path)
        visited.add(path)   
        files = (i for i in os.listdir(path) if i.endswith(".phy") or i.endswith('.phylip'))
        for file_name in files:
            if file_name not in visited_files:
                visited_files.add(file_name)
                phy_file_name = "{}/{}".format(path, file_name)
                with open("{}/{}".format(output_path, file_name), "w+") as new_handle:
                    with open(phy_file_name) as handle:
                        for line in handle:
                            
                            try:
                                split_line = line.split(" ")
                                if split_line[1] == "41":
                                    int_split = int(split_line[1])
                                    int_split -= len(remove_set)
                                    split_line[1] = " {} ".format(int_split)
                                    line = "".join((i for i in split_line))
                            except IndexError:
                                pass

                            if line.split(" ")[0] in remove_set:
                                continue
                            else:
                                new_handle.write(line)

        for root, subdirs, files in os.walk(path):
            for subdir in subdirs:
                return remove_ignored(subdir, output_path, visited)
    else:
        return True


def run_raxml(output_path):
    files = (("".join(file_name.rsplit(".")[:-1]), file_name) for file_name in os.listdir(output_path) if file_name.endswith(".phy") or file_name.endswith(".phylip")) 
    out_files = []
    good = set()
    commands = []

    for new_dir, new_fasta in files:
        output = "{}".format(new_dir)
        sp.call("mkdir -p {output_path}/{new_dir}".format(**locals()), shell=True)
        sp.call(("cp {}/{} ".format(output_path, new_fasta) + "{}/{}/{}".format(output_path, new_dir, new_fasta)), shell=True)
        random_number = 12345
        cmd = ['cd {output_path}/{new_dir} && raxmlHPC-AVX -f a  -m GTRGAMMA -x 12345 -# 500 -s {new_fasta} -n {output} -p {random_number} > RXML.{output}'.strip().format(**locals())]
        commands.append(cmd)
    return commands


def make_output_path(output_path, new_output="", count=0):
    if not os.path.isdir(output_path) and count == 0:
        os.makedirs(output_path)
        return output_path
    elif not os.path.isdir(new_output) and count > 0:
        os.makedirs(new_output)
        return new_output
    else:
        count += 1
        new_output = "{}_{}".format(output_path, count)
        return make_output_path(output_path, new_output, count)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-input_path', type=str, help="full_input_path", required=True)
    parser.add_argument('-output_path', type=str, help="full_output_path", required=True)
    parser.add_argument('-ignore', action='append', 
                                   dest='ignored',
                                   default=[],
                                   help='Ignore header add repeated values to a list',
                       )
    parser.add_argument('-c', type=int, default=1, dest='cores',help="number of cores you would like to use. default 1.")
    args = parser.parse_args()

    visited = set()
    visited_files = set()
    remove_set = set(args.ignored)

    output_path = make_output_path(args.output_path)
    status = remove_ignored(args.input_path, output_path)
    commands = run_raxml(output_path)

    threadList=[]
    for i in range(args.cores): #NUMBER CPU
        threadList.append("Thread-{}".format(i+1))

    queueLock = threading.Lock()
    workQueue = Queue.Queue(len(commands))
    threads = []
    threadID = 1

    exitFlag = 0

    # Create new threads
    for tName in threadList:
        thread = myThread(threadID, tName, workQueue)
        thread.start()
        threads.append(thread)
        threadID += 1

    # Fill the queue
    queueLock.acquire()
    for command in commands:
        workQueue.put(command)
    queueLock.release()


    # Wait for queue to empty
    while not workQueue.empty():
        pass
    # Notify threads it's time to exit
    exitFlag = 1

    # Wait for all threads to complete
    for thread in threads:
        thread.join()
    print("Exiting Main Thread")
    exitFlag = 0

    print("\nYOUR FILES HAVE BEEN CREATED AT \n\n{}\n".format(os.path.abspath(output_path)))
