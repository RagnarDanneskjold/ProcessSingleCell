import re
import operator

#class IntervalNode:
#    def __init__(self,int_start,int_end, int_index,sub_max):
#        self.leftChild = None
#        self.rightChild = None
#        self.start = int_start
#        self.end = int_end
#        self.index = int_index
#        self.subtreeMax = sub_max
#
#class IntervalTree:
#    def __init__(self):
#        self.root = None
#        tree_max = 0
#
#    def addNode(self, int_start, int_end, int_index):
#        if self.root == None:
#            self.root = IntervalNode(int_start, int_end, int_index,int_end)
#        else:
#            tree_max = self.add(int_start, int_end, int_index, self.root)
#
#    def compareNodes(self, node_start, node_end, int_start, int_end):
#        if int_start < node_start:
#            return -1 # left
#        elif int_start > node_start:
#            return 1 # right
#        elif int_start == node_start:
#            if int_end <= node_end:
#                return -1
#            else:
#                return 1
#
#    def add(self, int_start, int_end, int_index, node):
#        curr_node = node
#        print(curr_node.start, curr_node.end)
#        curr_max = int_end
#        while(1):
#            print(curr_node.start, curr_node.end)
#            if curr_max > curr_node.subtreeMax:
#                curr_node.subtreeMax = curr_max
#
#            if self.compareNodes(node.start, node.end, int_start, int_end) == -1:
#                if node.leftChild == None:
#                    print("adding left")
#                    node.leftChild = IntervalNode(int_start, int_end, int_index, curr_max)
#                    break
#                else:
#                    curr_node = node.leftChild
#            else:
#                if node.rightChild == None:
#                    print("adding right")
#                    node.rightChild = IntervalNode(int_start, int_end, int_index, curr_max)
#                    break
#                else:
#                    curr_node = node.rightChild
#        
#
#
#
#    def addRECURSIVE(self, int_start, int_end, int_index, node):
#        curr_max = int_end
##        print("node start:",node.start, "node.end:", node.end, "int start:", int_start, "int end:", int_end)
#        if self.compareNodes(node.start, node.end, int_start, int_end) == -1:
#            if node.leftChild == None:
#                node.leftChild = IntervalNode(int_start, int_end, int_index, curr_max)
#            else:
#                curr_max = self.add(int_start, int_end, int_index, node.leftChild)
#        else:
#            if node.rightChild == None:
#                node.rightChild = IntervalNode(int_start, int_end, int_index, curr_max)
#            else:
#                curr_max = self.add(int_start, int_end, int_index, node.rightChild)
#
#        if curr_max > node.subtreeMax:
#            node.subtreeMax = curr_max
#
#        return node.subtreeMax
#
#            
#    def nodeOverlap(self, node, int_start, int_end):
#        if node.start < int_end and node.end > int_start:
#            return True
#        else:
#            return False
#
#    def findNode(self, int_start, int_end):
#        if self.root == None:
#            return []
#        else:
#            return find(int_start, int_end, self.root, [])
#
#    def find(int_start, int_end, node, overlap_list):
#        curr_list = list(overlap_list)
#        if self.nodeOverlap(node, int_start, int_end):
#            curr_list = [overlap_list,node.index]
#
#        if (node.leftChild != None and node.subtreeMax > int_start): 
#            curr_list = self.find(int_start, int_end, node.leftChild, curr_list)
#
#        if node.rightChild != None:
#            curr_list = self.find(int_start, int_end, node.rightChild, curr_list)
#
#        return curr_list

class MyNode:
    def __init__(self, interval, gene_index, submin, submax):
        self.start = interval['start']
        self.end = interval['end']
        self.subtreeMax = submax
        self.subtreeMin = submin
        self.index = gene_index
        self.left = -1
        self.right = -1

class MyTree:
    def __init__(self):
        self.nodelist = []

    def comp(self, curr_node, interval):
        if interval['start'] < self.nodelist[curr_node].start:
            return -1
        elif interval['start'] > self.nodelist[curr_node].start:
            return 1
        else:
            if interval['end'] <= self.nodelist[curr_node].end:
                return -1
            else:
                return 1

    def addNode(self, interval, index):
        if self.nodelist == []:
            self.nodelist.append(MyNode(interval, index, interval['start'], interval['end']))
        else:
            curr_node = 0
            curr_max = self.nodelist[0].subtreeMax
            curr_min = self.nodelist[0].subtreeMin

            while(1):
                if curr_max > self.nodelist[curr_node].subtreeMax:
                    self.nodelist[curr_node].subtreeMax = curr_max

                if curr_min < self.nodelist[curr_node].subtreeMin:
                    self.nodelist[curr_node].subtreeMin = curr_min

                if self.comp(curr_node, interval) == -1:
                    if self.nodelist[curr_node].left == -1:
                        self.nodelist[curr_node].left = len(self.nodelist)
                        self.nodelist.append(MyNode(interval, index, curr_min, curr_max))
                        break
                    else:
                        curr_node = self.nodelist[curr_node].left
                else:
                    if self.nodelist[curr_node].right == -1:
                        self.nodelist[curr_node].right = len(self.nodelist)
                        self.nodelist.append(MyNode(interval, index, curr_min, curr_max))
                        break
                    else:
                        curr_node = self.nodelist[curr_node].right

    def nodeOverlap(self, interval, curr_node):
        if (self.nodelist[curr_node].start < interval['end'] and
                self.nodelist[curr_node].end > interval['start']):
            return True
        return False

    def findNode(self, interval):
#        print("nodelist:",len(self.nodelist))
        if self.nodelist == []:
            return []
        else:
            curr_node = 0
            overlap_list = []
            node_stack = []
            stack_curr_index = -1
            done = False
            while not done:
                if curr_node != -1:
                    node_stack.append(curr_node)
                    if (self.nodelist[curr_node].left != -1 and
                            self.nodelist[self.nodelist[curr_node].left].subtreeMax > interval['start']):
                        curr_node = self.nodelist[curr_node].left
                    else:
                        curr_node = -1
#                    print("left:", curr_node)
                else:
                    if len(node_stack) > 0:
                        curr_node = node_stack.pop()
#                        print(curr_node)
                        
                        if self.nodeOverlap(interval, curr_node):
                            overlap_list.append(self.nodelist[curr_node].index)

                        if (self.nodelist[curr_node].right != -1 and
                                self.nodelist[self.nodelist[curr_node].right].subtreeMin < interval['end']):
                            curr_node = self.nodelist[curr_node].right
                        else:
                            curr_node = -1
                    else:
                        done = True

            return overlap_list

                
def parseValidLine(line):
    line_ar = line.rstrip('\n\r').split('\t')

    gene_of_interest = False
#    if re.search(r'gene_id \"ENSMUSG00000028894\.18',line):
#        print ("ON GENE OF INTEREST")
#        gene_of_interest = True

    if len(line_ar) < 9:
#        print("line < 9")
        return {}

    id_string = line_ar[8]

    gene_id = re.search(r'gene_id \"(.+?)\";', id_string).group(1)

#    if gene_of_interest:
#        print("GeneID:",gene_id)

    result = {
            'chrom': line_ar[0],
            'feature': line_ar[2],
            'start': int(line_ar[3]),
            'end': int(line_ar[4]),
            'strand': line_ar[6],
            'gene_id': gene_id
        }

#    if gene_of_interest:
#        print("feature:",result['feature'])
    if result['feature'] != "gene" and result['feature'] != "exon":
 #       if gene_of_interest:
  #          print("BAD FEATURE:{}")
        return {}

#    if gene_of_interest:
#        print ("returning result")
#        print(result)

    return result

#def sortGTFStructure(data):
#    for chrom in data:
#        for strand in data[chrom]:
#            data[chrom][strand].sort(key=operator.itemgetter('start', 'end'))
#            for i in range(len(data[chrom][strand])):
#                data[chrom][strand][i]['exons'].sort(key=operator.itemgetter('start', 'end'))

#def checkForGene(data):
#    for chrom in data:
#        for strand in data[chrom]:
#            for i in range(len(data[chrom][strand])):
#                if data[chrom][strand][i]['gene_id'] == 'ENSMUSG00000028894.18':
#                    print("found gene of interest at index",i,"chrom",chrom,"strand:",strand)
#                    print ("i-1:",data[chrom][strand][i-1])
#                    print ("i:",data[chrom][strand][i])
#                    print ("i+1:",data[chrom][strand][i+1])

def parseGTFFile (gtf_fp):
    parsedData = dict()

    count = 0
    for line in gtf_fp:
        count+=1
        if line.startswith('#'):
            continue

        fields = parseValidLine(line)

        if not fields:
            continue
        
        if fields['chrom'] not in parsedData:
            parsedData[fields['chrom']] = {
                '+': {
                    'tree' : MyTree(),
                    'genes' : list()
                    },
                '-': {
                    'tree' : MyTree(),
                    'genes' : list()
                    }
                }

        if fields['feature'] == 'gene':
            curr_index = len(parsedData[fields['chrom']][fields['strand']]['genes'])
#            if fields['gene_id'] == 'ENSMUSG00000028894.18':
#                print("inserting gene of interest into structure")

            parsedData[fields['chrom']][fields['strand']]['tree'].addNode(fields, curr_index)

            parsedData[fields['chrom']][fields['strand']]['genes'].append({
                'gene_id':fields['gene_id'],
                'start':fields['start'],
                'end':fields['end'],
                'exons':list(),
                'reads':list()
                })

        else: # exon
            curr_index = len(parsedData[fields['chrom']][fields['strand']]['genes']) - 1
            parsedData[fields['chrom']][fields['strand']]['genes'][curr_index]['exons'].append({
                'start':fields['start'],
                'end':fields['end']
                })
      
    return parsedData
