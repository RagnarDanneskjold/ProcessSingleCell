class MyNode:
    def __init__(self, interval, gene_index, submin, submax):
        self.start = interval['start']
        self.end = interval['end']
        self.subtreeMax = submax
        self.subtreeMin = submin
        self.index = gene_index
        self.left = -1
        self.right = -1

    def getIndex(self):
        return self.index

class MyTree:
    def __init__(self):
        self.nodelist = []
        self.root = -1

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

#    def addNode(self, interval, index):
    def addNode(self, interval):
        if self.nodelist == []:
#            self.nodelist.append(MyNode(interval, index, interval['start'], interval['end']))
            self.root = 0
            self.nodelist.append(MyNode(interval, self.root, interval['start'], interval['end']))
        else:
            curr_node = self.root
            curr_max = self.nodelist[self.root].subtreeMax
            curr_min = self.nodelist[self.root].subtreeMin

            while(1):

                if self.comp(curr_node, interval) == -1:
                    if interval['start'] < self.nodelist[curr_node].subtreeMin:
                        self.nodelist[curr_node].subtreeMin = interval['start']

                    if self.nodelist[curr_node].left == -1:
                        my_index = len(self.nodelist)
                        self.nodelist[curr_node].left = my_index
                        self.nodelist.append(MyNode(interval, my_index, curr_min, curr_max))
                        break
                    else:
                        curr_node = self.nodelist[curr_node].left
                else:
                    if interval['end'] > self.nodelist[curr_node].subtreeMax:
                        self.nodelist[curr_node].subtreeMax = interval['end']
                    if self.nodelist[curr_node].right == -1:
                        my_index = len(self.nodelist)
                        self.nodelist[curr_node].right = my_index
                        self.nodelist.append(MyNode(interval, my_index, curr_min, curr_max))
                        break
                    else:
                        curr_node = self.nodelist[curr_node].right

    def recursive_rebuild(self, sort_nodes, start_index, end_index):
        if start_index > end_index:
            return []

        mid_index = (start_index + end_index) // 2
        curr_root = sort_nodes[mid_index]

        left_values = self.recursive_rebuild(sort_nodes, start_index, mid_index - 1)
        if left_values != []:
            self.nodelist[curr_root].left = left_values[0] 

            if left_values[1] < self.nodelist[curr_root].start:
                self.nodelist[curr_root].subtreeMin = left_values[1] 
            else:
                self.nodelist[curr_root].subtreeMin = self.nodelist[curr_root].start 

        else:
            self.nodelist[curr_root].left = -1 
            self.nodelist[curr_root].subtreeMin = self.nodelist[curr_root].start

        right_values = self.recursive_rebuild(sort_nodes, mid_index + 1, end_index)
        if right_values != []:
            self.nodelist[curr_root].right = right_values[0]

            if right_values[2] > self.nodelist[curr_root].end:
                self.nodelist[curr_root].subtreeMax = right_values[2]
            else:
                self.nodelist[curr_root].subtreeMax = self.nodelist[curr_root].end
        else:
            self.nodelist[curr_root].right = -1
            self.nodelist[curr_root].subtreeMax = self.nodelist[curr_root].end


        return [curr_root, self.nodelist[curr_root].subtreeMin, self.nodelist[curr_root].subtreeMax]

    def balance(self):
        # perform in-order traversal
        sort_nodes = []
        if self.nodelist == []:
            return 
        else:
            curr_node = self.root
            node_stack = []
            done = False
            while not done:
                if curr_node != -1:
                    node_stack.append(curr_node)
                    curr_node = self.nodelist[curr_node].left
                else:
                    if len(node_stack) > 0:
                        curr_node = node_stack.pop()
            
                        sort_nodes.append(curr_node)

                        curr_node = self.nodelist[curr_node].right
                    else:
                        done = True

            
        # now rebuild tree
        root_values = self.recursive_rebuild(sort_nodes, 0,len(sort_nodes)-1)
        self.root = root_values[0]

#        print("ROOT:",self.root)
#        for i in range(len(self.nodelist)):
#            print("BALANCE: root:",i,"left:",self.nodelist[i].left,"right:",self.nodelist[i].right, "my_min:",self.nodelist[i].subtreeMin, "my_max:",self.nodelist[i].subtreeMax, "start:",self.nodelist[i].start, "end:",self.nodelist[i].end)
            


    def nodeOverlap(self, interval, curr_node):
        if (self.nodelist[curr_node].start < interval['end'] and
                self.nodelist[curr_node].end > interval['start']):
            return True
        return False

class ExonTree(MyTree):
    def __init__(self):
        MyTree.__init__(self)

    def findNodeBool(self, interval):
        if self.nodelist == []:
            return False
        else:
            curr_node = self.root
            node_stack = []
            done = False

            while not done:
                if curr_node != -1:
                    node_stack.append(curr_node)

                    if self.nodelist[self.nodelist[curr_node].left].subtreeMax > interval['start']:
                        curr_node = self.nodelist[curr_node].left
                    else:
                        curr_node = -1
                else:
                    if len(node_stack) > 0:
                        curr_node = node_stack.pop()

                        if self.nodeOverlap(interval, curr_node):
                            return True

                        if (self.nodelist[curr_node].right != -1 and
                                self.nodelist[self.nodelist[curr_node].right].subtreeMin < interval['end']):
                            curr_node = self.nodelist[curr_node].right
                        else:
                            curr_node = -1
                    else:
                        done = True

            return False

class GeneNode(MyNode):
    def __init__(self,interval,geneid,gene_index,submin,submax):
        MyNode.__init__(self,interval,gene_index,submin,submax)
        self.gene_id = geneid
        self.exons = ExonTree()
        self.reads = 0

    def getID(self):
        return self.gene_id

    def getReads(self):
        return self.reads

    def checkExons(self, interval):
        return self.exons.findNodeBool(interval)

    def writeOut(self, out_fp, mychrom, mystrand):
        outstr = mychrom + "\t" + str(self.start) + "\t" + str(self.end) + "\t" + mystrand + "\t" + self.gene_id + "\t" + str(self.reads) + "\n"
        out_fp.write(outstr)

class GeneTree(MyTree):
    def __init__(self, chrom, strand):
        MyTree.__init__(self)
        self.chrom = chrom
        self.strand = strand

    def addExon(self, fields, gene_index):
        self.nodelist[gene_index].exons.addNode(fields)


    def strictOverlap(self,interval,index):
        if (interval['start'] >= self.nodelist[index].start and
                interval['end'] <= self.nodelist[index].end):
            return True
        return False

#    def linearSearch(self,interval):
#        for i in range(len(self.nodelist)):
#            if self.nodeOverlap(interval,i):
#                print("LINEAR_SEARCH_OVERLAP: [",i,"]:", self.nodelist[i].start, self.nodelist[i].end, interval)

    def overlapInterval(self,interval, TMP_PREV_READS):
        read = TMP_PREV_READS[0]
        print("Checking Overlap:Reads:\t", TMP_PREV_READS[0].query_name, "\t", TMP_PREV_READS[1].query_name)
        indicies = self.findNode(interval)

#        if TMP_PREV_READS[0].query_name == "HWI-ST999:184:C44V8ACXX:7:1101:1362:54252":
#            self.linearSearch(interval)
#            print("linear search done")

        if (indicies == []): #or 
                #len(indicies) > 1 or
#                not self.strictOverlap(interval, indicies[0])
#                ):
#            if indicies == []:
            print("not_overlap! interval:",interval)
#            elif len(indicies) > 1:
#                for i in indicies:
#                    print("Gene i:",self.nodelist[i].gene_id, "start:", self.nodelist[i].start, "end:", self.nodelist[i].end)
#                print("multiple_overlap!")
#            else:
#                print("not_strict_overlap!")
            print("UNASSIGNED:",read.query_name)
            return False

        count_exon_overlaps = 0
        true_genes = []
        if len(indicies) > 1:
            for i in indicies:
                print("multiple indicies i:",i)
                print("Gene", i,":",self.nodelist[i].gene_id, "start:", self.nodelist[i].start, "end:", self.nodelist[i].end)
                if (self.nodelist[i].checkExons(interval)):# and 
#                        self.strictOverlap(interval,i)):
                    count_exon_overlaps += 1
                    true_genes.append(i)
            print("multiple_overlap!:exon overlaps:",count_exon_overlaps)

            if count_exon_overlaps > 1:
                print("exon overlaps > 1:")
                for i in true_genes:
                    print("Exon overlap Gene i:",self.nodelist[i].gene_id, "start:", self.nodelist[i].start, "end:", self.nodelist[i].end)
                print("UNASSIGNED:",read.query_name)
                return False
            elif count_exon_overlaps == 0:
                print("UNASSIGNED:",read.query_name)
                print("exon overlaps == 0")
                return False
            elif count_exon_overlaps == 1:
                print("exon overlaps == 1! SUCCESS Adding read")
                self.nodelist[true_genes[0]].reads += 1
                return True

        
        if (self.nodelist[indicies[0]].checkExons(interval)):# and 
#                self.strictOverlap(interval,indicies[0])):
            print("SUCCESS_read_added!")

#            if self.nodelist[indicies[0]].getID() == "ENSMUSG00000062794.8":
#                print("ID:[",self.nodelist[indicies[0]].getID(),"] start:",self.nodelist[indicies[0]].start, "end:", self.nodelist[indicies[0]].end)
#                print("Reads:\t", TMP_PREV_READS[0].query_name, "\t", TMP_PREV_READS[1].query_name)


            self.nodelist[indicies[0]].reads += 1
            return True

        if not self.nodelist[indicies[0]].checkExons(interval):
            print("exon check failed!")
#        if not self.strictOverlap(interval,indicies[0]):
#            print("strict overlap failed!")
        print("indicies len is:",len(indicies))
#       print("no_exon_overlap!")
        print("UNASSIGNED:",read.query_name)
        return False

    def writeTree(self, out_fp, mychrom, mystrand):
        for i in self.nodelist: # UNSORTED OUTPUT
            i.writeOut(out_fp, mychrom, mystrand)

    def findNode(self, interval):
        if self.nodelist == []:
            return []
        else:
            curr_node = self.root
            overlap_list = []
            node_stack = []
            done = False
            while not done:
#                print("curr node:",curr_node)
                if curr_node != -1:
#                    print("NODE != -1: stack pre append:",node_stack)
                    node_stack.append(curr_node)
                    if (self.nodelist[curr_node].left != -1 and
                            self.nodelist[self.nodelist[curr_node].left].subtreeMax > interval['start']):
                        curr_node = self.nodelist[curr_node].left
                    else:
                        curr_node = -1
                else:
#                    print("NODE == -1: stack pre pop:",node_stack)
                    if len(node_stack) > 0:
                        curr_node = node_stack.pop()
                        
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


    def addNode(self, interval, gene_id):
        if self.nodelist == []:
            self.root = 0
            self.nodelist.append(GeneNode(interval, gene_id, self.root, interval['start'], interval['end']))
            return self.root
        else:
            curr_node = self.root
            curr_max = self.nodelist[self.root].subtreeMax
            curr_min = self.nodelist[self.root].subtreeMin

            while(1):

                if self.comp(curr_node, interval) == -1:
                    if interval['start'] < self.nodelist[curr_node].subtreeMin:
                        self.nodelist[curr_node].subtreeMin = interval['start']

                    if self.nodelist[curr_node].left == -1:
                        my_index = len(self.nodelist)
                        self.nodelist[curr_node].left = my_index
                        self.nodelist.append(GeneNode(interval, gene_id, my_index, curr_min, curr_max))
                        return my_index
                    else:
                        curr_node = self.nodelist[curr_node].left
                else:
                    if interval['end'] > self.nodelist[curr_node].subtreeMax:
                        self.nodelist[curr_node].subtreeMax = interval['end']
                    if self.nodelist[curr_node].right == -1:
                        my_index = len(self.nodelist)
                        self.nodelist[curr_node].right = my_index
                        self.nodelist.append(GeneNode(interval, gene_id, my_index, curr_min, curr_max))
                        return my_index
                    else:
                        curr_node = self.nodelist[curr_node].right

    def balanceAll(self):
        self.balance()

        for gene in self.nodelist:
            gene.exons.balance()
