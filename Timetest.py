# Aparna Rajpurkar
# Timetest script to compare runtimes of HTseq, FeatureCounts, and FastCount
# uses a list of chromosome-split BAM files as input

# imports
import sys
import subprocess
import re

## CONSTANT ##
REPS = 3

## FUNCTIONS ##

def parseTime(time_output):
    """Parses output of UNIX time program to seconds"""
    # use regex to get time in minutes and seconds
    get_time = re.compile(r'system\s(.+)elapsed')
    match = get_time.search(time_output).groups(1)
    # split minutes and seconds by colon
    time = match[0].split(':')
    # calculate seconds
    time_in_secs = float(time[0])*60+float(time[1])
    
    # return time in seconds
    return time_in_secs

def runFastCount(bams,cmd):
    """Run the FastCount program"""
    # set tmpfile name and open it
    tmpfile = "tmp"
    tmp_fp = open(tmpfile,"w")
    # print bam list to tmpfile as required by FastCount
    tmp_fp.write("\n".join(bams))
    # close tmpfile
    tmp_fp.close()
    # run FastCount and get time output
    time_output = subprocess.getoutput(cmd + tmpfile)
    # parse it and return time in seconds
    return parseTime(time_output)

def runFeatureCounts(bams,cmd):
    """Run FeatureCounts program on input bams"""
    # run FeatureCounts program
    time_output = subprocess.getoutput(cmd + " ".join(bams))
    # return time in seconds
    return parseTime(time_output)

def runHTSeq(bams,cmd1,cmd2,outfile_base):
    """run HTSeq program on input bams"""

    # HTSeq cannot handle multiple BAM inputs, so we run it individually 
    # for each, and add up the time
    # initialize variables
    count = 0
    curr_time = 0

    # iterate over each bamfile
    for bam in bams:
        # construct command string and tmp output file
        cmd = cmd1 + bam + cmd2 + outfile + str(count) + ".htseq"
        # run command
        time_output = subprocess.getoutput(cmd)
        # increment count
        count += 1
        # get time in seconds and add to current time
        curr_time += parseTime(time_output)

    # return time in seconds
    return curr_time

# parse command line input. 
if len(sys.argv) < 4:
    # check if args are correct number
    # print usage and exit if not
    print("usage: <gtffile> <outfile basename counts> <full_bam> <bamfile_list> <time outfile>", file=sys.stderr)
    sys.exit(1)

# parse sys.argv to useful variables
gtffile = sys.argv[1]
outfile = sys.argv[2]
bamfile = sys.argv[3]
bamfile_list_chroms = sys.argv[4]
time_outfile = sys.argv[4]

# open list of chromosomes
bam_fp = open(bamfile_list_chroms, "r")
bams = []

# read in list of chromosomes, append to bams list
for line in bam_fp:
    parsed_line = line.rstrip('\n\r')
    bams.append(parsed_line)

# close file
bam_fp.close()

# construct commands for each program
cmd_FastCount = "time python3 FastCount.py " + gtffile + " " + outfile + ".FastCount " 
cmd_FeatureCounts = "time featureCounts -a " + gtffile + " -o " + outfile + ".FeatureCount" + " -s 1 -p -B -C -R --primary " 
cmd_HTseq_1 = "time htseq-count -f bam -r name -s yes -t gene "
cmd_HTseq_2 = " " + gtffile + " > " 

# initialize time dictionary
# that keeps track of times for each test
times = {
        'fastcount': list(),
        'featurec': list(),
        'htseq': list(),
        'numfiles': list()
        }

# initialize num_files, loop variable
num_files = 1
last_loop = False

# while number of files to test is less than total number of files and 
# we want to keep looping
while num_files <= len(bams) or last_loop:
    # get the bam names corresponding to this test
    curr_bams = bams[:num_files]

    # if last loop, we want to run the test on the FULL bamfile, not split
    if last_loop:
        num_files = "full"
        curr_bams = bamfile

    # initialize times to 0
    time_fa = 0
    time_fc = 0
    time_ht = 0

    # repeat each test 3 times
    for i in range(REPS):
        # run each program and capture the time
        time_fa += runFastCount(curr_bams,cmd_FastCount)
        time_fc += runFeatureCounts(curr_bams,cmd_FeatureCounts)
        time_ht += runHTSeq(curr_bams, cmd_HTseq_1, cmd_HTseq_2, outfile)

    # calculate average time for each program and add to dictionary
    times['fastcount'].append(str(time_fa/REPS))
    times['featurec'].append(str(time_fc/REPS))
    times['htseq'].append(str(time_ht/REPS))
    # add number of files for this test to the dictionary
    times['numfiles'].append(str(num_files))

    # construct a header line and all output lines for this test
    header = "Program\t" + "\t".join(times['numfiles']) + "\n"
    fa = "FastCount\t" + "\t".join(times['fastcount']) + "\n"
    fc = "FeatureCounts\t" + "\t".join(times['featurec']) + "\n"
    ht = "HTSeq\t" + "\t".join(times['htseq']) + "\n"

    # output results for this test
    out_fp = open(time_outfile + "_NUMFILES_" + str(num_files),"w")

    out_fp.write(header)
    out_fp.write(fa)
    out_fp.write(fc)
    out_fp.write(ht)

    out_fp.close()

    if not last_loop:
        # if not last loop, we double the test
        num_files *= 2

        # if this is too many files, set flag to do the final test and exit
        if num_files >= len(bams):
            last_loop = True
    else:
        last_loop = True

