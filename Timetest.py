import sys
import subprocess
import re

def parseTime(time_output):
    get_time = re.compile(r'system\s(.+)elapsed')
    match = get_time.search(time_output).groups(1)
    time = match[0].split(':')
    time_in_secs = float(time[0])*60+float(time[1])

    return time_in_secs

def runFastCount(bams,cmd):
    print("Fast Count",file=sys.stderr)
    tmpfile = "tmp"
    tmp_fp = open(tmpfile,"w")
    tmp_fp.write("\n".join(bams))
    tmp_fp.close()
    time_output = subprocess.getoutput(cmd + tmpfile)
    return parseTime(time_output)

def runFeatureCounts(bams,cmd):
    print("Feature Count",file=sys.stderr)
    time_output = subprocess.getoutput(cmd + " ".join(bams))
    return parseTime(time_output)

def runHTSeq(bams,cmd1,cmd2,outfile_base):
    print("HTSeq",file=sys.stderr)
    count = 0
    curr_time = 0

    for bam in bams:
        cmd = cmd1 + bam + cmd2 + outfile + str(count) + ".htseq"
        time_output = subprocess.getoutput(cmd)
        count += 1
        curr_time += parseTime(time_output)

    return curr_time

if len(sys.argv) < 4:
    print("usage: <gtffile> <outfile basename counts> <full_bam> <bamfile_list> <time outfile>", file=sys.stderr)
    sys.exit(1)

gtffile = sys.argv[1]
outfile = sys.argv[2]
bamfile = sys.argv[3]
bamfile_list_chroms = sys.argv[4]
time_outfile = sys.argv[4]

bam_fp = open(bamfile_list_chroms, "r")
bams = []

for line in bam_fp:
    parsed_line = line.rstrip('\n\r')
    bams.append(parsed_line)

bam_fp.close()

cmd_FastCount = "time python3 FastCount.py " + gtffile + " " + outfile + ".FastCount " 
cmd_FeatureCounts = "time featureCounts -a " + gtffile + " -o " + outfile + ".FeatureCount" + " -s 1 -p -B -C -R --primary " 
cmd_HTseq_1 = "time htseq-count -f bam -r name -s yes -t gene "
cmd_HTseq_2 = " " + gtffile + " > " 

#time = runHTSeq(bams,cmd_HTseq_1,cmd_HTseq_2,outfile)
#time = runFeatureCounts(bams,cmd_FeatureCounts)
#print("TIME OUTPUT:[",time,"]")

times = {
        'fastcount': list(),
        'featurec': list(),
        'htseq': list(),
        'numfiles': list()
        }

num_files = 1
REPS = 3
last_loop = False

while num_files <= len(bams) or last_loop:
    print("Numfiles is:",num_files,file=sys.stderr)
    curr_bams = bams[:num_files]

    if last_loop:
        num_files = "full"
        curr_bams = bamfile

    time_fa = 0
    time_fc = 0
    time_ht = 0

    for i in range(REPS):
        time_fa += runFastCount(curr_bams,cmd_FastCount)
        time_fc += runFeatureCounts(curr_bams,cmd_FeatureCounts)
        time_ht += runHTSeq(curr_bams, cmd_HTseq_1, cmd_HTseq_2, outfile)

    times['fastcount'].append(str(time_fa/REPS))
    times['featurec'].append(str(time_fc/REPS))
    times['htseq'].append(str(time_ht/REPS))
    times['numfiles'].append(str(num_files))

    header = "Program\t" + "\t".join(times['numfiles']) + "\n"
    fa = "FastCount\t" + "\t".join(times['fastcount']) + "\n"
    fc = "FeatureCounts\t" + "\t".join(times['featurec']) + "\n"
    ht = "HTSeq\t" + "\t".join(times['htseq']) + "\n"

    out_fp = open(time_outfile + "_NUMFILES_" + str(num_files),"w")

    out_fp.write(header)
    out_fp.write(fa)
    out_fp.write(fc)
    out_fp.write(ht)

    out_fp.close()

    if not last_loop:
        num_files *= 2

        if num_files >= len(bams):
            last_loop = True
    else:
        last_loop = True

