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

    def addNode(self, interval):
        if self.nodelist == []:
            self.root = 0
            self.nodelist.append(
                    MyNode(interval, self.root, interval['start'], interval['end'])
                    )
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
                        self.nodelist.append(
                                MyNode(interval, my_index, curr_min, curr_max)
                                )
                        break
                    else:
                        curr_node = self.nodelist[curr_node].left
                else:
                    if interval['end'] > self.nodelist[curr_node].subtreeMax:
                        self.nodelist[curr_node].subtreeMax = interval['end']
                    if self.nodelist[curr_node].right == -1:
                        my_index = len(self.nodelist)
                        self.nodelist[curr_node].right = my_index
                        self.nodelist.append(
                                MyNode(interval, my_index, curr_min, curr_max)
                                )
                        break
                    else:
                        curr_node = self.nodelist[curr_node].right

    def rebuild(self, sort_nodes, start_index, end_index):
        if start_index > end_index:
            return []

        mid_index = (start_index + end_index) // 2
        curr_root = sort_nodes[mid_index]

        left_values = self.rebuild(sort_nodes, start_index, mid_index - 1)
        if left_values != []:
            self.nodelist[curr_root].left = left_values[0] 

            if left_values[1] < self.nodelist[curr_root].start:
                self.nodelist[curr_root].subtreeMin = left_values[1] 
            else:
                self.nodelist[curr_root].subtreeMin = self.nodelist[curr_root].start 

        else:
            self.nodelist[curr_root].left = -1 
            self.nodelist[curr_root].subtreeMin = self.nodelist[curr_root].start

        right_values = self.rebuild(sort_nodes, mid_index + 1, end_index)
        if right_values != []:
            self.nodelist[curr_root].right = right_values[0]

            if right_values[2] > self.nodelist[curr_root].end:
                self.nodelist[curr_root].subtreeMax = right_values[2]
            else:
                self.nodelist[curr_root].subtreeMax = self.nodelist[curr_root].end
        else:
            self.nodelist[curr_root].right = -1
            self.nodelist[curr_root].subtreeMax = self.nodelist[curr_root].end


        min_val = self.nodelist[curr_root].subtreeMin
        max_val = self.nodelist[curr_root].subtreeMax
        return [curr_root, min_val, max_val]

    def inOrderTraverse(self):
        sort_nodes = []
        if self.nodelist == []:
            return [] 
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

        return sort_nodes    

    def balance(self):
        # perform in-order traversal
        sort_nodes = self.inOrderTraverse()
        if sort_nodes == []:
            return
        # now rebuild tree
        root_values = self.rebuild(sort_nodes, 0,len(sort_nodes)-1)
        self.root = root_values[0]
            
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

                    left = self.nodelist[curr_node].left
                    if self.nodelist[left].subtreeMax > interval['start']:
                        curr_node = self.nodelist[curr_node].left
                    else:
                        curr_node = -1
                else:
                    if len(node_stack) > 0:
                        curr_node = node_stack.pop()

                        if self.nodeOverlap(interval, curr_node):
                            return True

                        right = self.nodelist[curr_node].right
                        if (right != -1 and
                                self.nodelist[right].subtreeMin < interval['end']):
                            curr_node = self.nodelist[curr_node].right
                        else:
                            curr_node = -1
                    else:
                        done = True

            return False

class GeneNode(MyNode):
    def __init__(self,interval,geneid,gene_index,submin,submax,num_of_bams):
        MyNode.__init__(self,interval,gene_index,submin,submax)
        self.gene_id = geneid
        self.exons = ExonTree()
        self.reads = [0] * num_of_bams

    def addReads(self, bam_num):
        self.reads[bam_num] += 1

    def checkExons(self, interval):
        return self.exons.findNodeBool(interval)

    def writeOut(self, out_fp, mychrom, mystrand):
        outstr = mychrom + "\t" + str(self.start) + "\t" + str(self.end) +\
                "\t" + mystrand + "\t" + self.gene_id 
        
        for i in self.reads:
            outstr += "\t" + str(i)

        outstr += "\n"

        out_fp.write(outstr)

class GeneTree(MyTree):
    def __init__(self, chrom, strand):
        MyTree.__init__(self)
        self.chrom = chrom
        self.strand = strand

    def addExon(self, fields, gene_index):
        self.nodelist[gene_index].exons.addNode(fields)

    def overlapInterval(self,interval, bam_num):
        indicies = self.findNode(interval)

        if indicies == []: 
            return False

        count_exon_overlaps = 0
        true_gene = -1

        if len(indicies) > 1:
            for i in indicies:
                if self.nodelist[i].checkExons(interval): 
                    count_exon_overlaps += 1
                    true_gene = i

            if count_exon_overlaps > 1 or count_exon_overlaps == 0:
                return False
            elif count_exon_overlaps == 1:
                self.nodelist[true_gene].addReads(bam_num)
                return True

        
        if self.nodelist[indicies[0]].checkExons(interval): 
            self.nodelist[indicies[0]].addReads(bam_num)
            return True

        return False

    def writeTree(self, out_fp):
        sort_nodes = self.inOrderTraverse()

        for i in sort_nodes: 
            self.nodelist[i].writeOut(out_fp, self.chrom, self.strand)

    def findNode(self, interval):
        if self.nodelist == []:
            return []
        else:
            curr_node = self.root
            overlap_list = []
            node_stack = []
            done = False

            while not done:
                if curr_node != -1:
                    node_stack.append(curr_node)
                    left = self.nodelist[curr_node].left
                    if (left != -1 and
                            self.nodelist[left].subtreeMax > interval['start']):
                        curr_node = self.nodelist[curr_node].left
                    else:
                        curr_node = -1
                else:
                    if len(node_stack) > 0:
                        curr_node = node_stack.pop()
                        
                        if self.nodeOverlap(interval, curr_node):
                            overlap_list.append(self.nodelist[curr_node].index)

                        right = self.nodelist[curr_node].right
                        if (right != -1 and
                                self.nodelist[right].subtreeMin < interval['end']):
                            curr_node = self.nodelist[curr_node].right
                        else:
                            curr_node = -1
                    else:
                        done = True

            return overlap_list


    def addNode(self, interval, gene_id,bam_num):
        if self.nodelist == []:
            self.root = 0
            self.nodelist.append(
                    GeneNode(interval, gene_id, self.root, 
                        interval['start'], interval['end'], bam_num)
                    )
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
                        self.nodelist.append(
                                GeneNode(interval, gene_id, my_index, 
                                    curr_min,curr_max, bam_num)
                                )
                        return my_index
                    else:
                        curr_node = self.nodelist[curr_node].left
                else:
                    if interval['end'] > self.nodelist[curr_node].subtreeMax:
                        self.nodelist[curr_node].subtreeMax = interval['end']
                    if self.nodelist[curr_node].right == -1:
                        my_index = len(self.nodelist)
                        self.nodelist[curr_node].right = my_index
                        self.nodelist.append(
                                GeneNode(interval, gene_id, my_index, 
                                    curr_min, curr_max,bam_num)
                                )
                        return my_index
                    else:
                        curr_node = self.nodelist[curr_node].right

    def balanceAll(self):
        self.balance()

        for gene in self.nodelist:
            gene.exons.balance()
