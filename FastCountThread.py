import sys
import pysam
import Bam, GTFparse, InputOutput, Reads, Tree
from threading import Thread

#class MyThread(threading.Thread):
#    def __init__(self, file_list, i_list, genes):
#        threading.Thread.__init__(self)
#
#        self.files = file_list
#        self.indicies = i_list
#        self.genes = genes
#
#    def run(self):
#        for index in range(len(self.files)):
#            print("Start BAM file [",self.indicies[index],"]", file=sys.stderr)
#            Bam.processBam(self.files[index],self.indicies[index],genes)
#            print("End BAM file [",self.indicies[index],"]", file=sys.stderr)
#
THREADING = True

files = InputOutput.getInput(sys.argv,THREADING)

if not files:
    print("usage: python3 FastCount.py <path/to/GTF> <path/to/outfile> <number of threads> <path/to/bamfile_1> <path/to/bamfile_2> <...>")
    sys.exit(1)

thread_num = files['thread']
bams_per_thread = int(len(files['bams']) / int(thread_num))
print("bams per thread",bams_per_thread)

print("Start GTF file", file=sys.stderr)
genes = GTFparse.parseGTFFile(files['gtffile'], len(files['bams']))
print("End GTF file", file=sys.stderr)




def threadBams(filelist,i_list,genes):
    for index in range(len(filelist)):
        print("Start BAM file [",i_list[index],"]", file=sys.stderr)
        Bam.processBam(filelist[index],i_list[index],genes)
        print("Start BAM file [",i_list[index],"]", file=sys.stderr)
    return

thread_list = [] * int(thread_num)
tmpfilelist = list()
i_list = list()

for i in range(len(files['bams'])):
    if i % bams_per_thread == 0:
        t = Thread(target=threadBams, args=(tmpfilelist,i_list,genes,))
#        thread_list[len(thread_list)-1].start()
        t.start()
        thread_list.append(t)
        tmpfilelist = []
        i_list = []

    i_list.append(i)
    tmpfilelist.append(files['bams'][i])
    
if i_list != []:
    t = Thread(target=threadBams, args=(tmpfilelist,i_list,genes,))
    t.start()
    thread_list.append(t)
    tmpfilelist = []
    i_list = []

for t in thread_list:
    t.join()

print("Start outputting", file=sys.stderr)
InputOutput.printOutput(files['outfile'],genes,files['bams'])
print("End outputting", file=sys.stderr)
