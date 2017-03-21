# Aparna Rajpurkar
# Tree class definitions for FastCount.py

## BASE CLASSES ##
class MyNode:
    """Node Base class. This holds information about the gene/exon object"""
    def __init__(self, interval, gene_index, submin, submax):
        """Initialization function: set information about this object\
                and initialize it's left and right children to NULL equivalent"""
        # node information
        self.start = interval['start']
        self.end = interval['end']
        self.subtreeMax = submax
        self.subtreeMin = submin
        self.index = gene_index

        # child node links
        self.left = -1
        self.right = -1

class MyTree:
    """Base Tree class, inherited in ExonTree and GeneTree classes. Allows\
            tree traversal, node comparisons, and balance operations"""

    def __init__(self):
        """ initialization function : initialize the tree's node list and set\
                root to NULL equivalent"""
        self.nodelist = []
        self.root = -1

    def comp(self, curr_node, interval):
        """compare an interval to a node, return boolean"""
        
        # first compare start sites, if equal compare end sites.
        # ties are broken with returning -1
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
        """Basic node adding function. Take an interval dictionary, find \
                correct place in tree to add it using an iterative operation,\
                add node"""

        if self.nodelist == []:
            # if we have no nodes in our Tree, add this one as root
            self.root = 0
            self.nodelist.append(
                    MyNode(interval, self.root, interval['start'], interval['end'])
                    )
        else:
            # initialize current node variables
            curr_node = self.root
            curr_max = self.nodelist[self.root].subtreeMax
            curr_min = self.nodelist[self.root].subtreeMin

            # while we have not broken out of the loop, search
            # the tree using an iterative method
            # because the trees get very large, cannot use recursive
            while(1):
                # check if our interval is less than curr_node
                if self.comp(curr_node, interval) == -1:
                    # check left node
                    # check if the current subtreeMin of the left child is
                    # less than our start
                    if interval['start'] < self.nodelist[curr_node].subtreeMin:
                        # if our start is smaller, it becomes the new subtreeMin
                        self.nodelist[curr_node].subtreeMin = interval['start']

                    # check if the left child is empty
                    if self.nodelist[curr_node].left == -1:
                        # if it's empty, add a new node: this interval
                        # as the left child
                        # get current index of nodelist where we will add this node
                        my_index = len(self.nodelist)
                        # set the curr_node left child link to our index
                        self.nodelist[curr_node].left = my_index
                        # append our interval as a new node
                        self.nodelist.append(
                                MyNode(interval, my_index, curr_min, curr_max)
                                )
                        # done! break out of loop
                        break
                    else:
                        # if left child is not empty, search the left subtree
                        # on the next loop
                        curr_node = self.nodelist[curr_node].left
                else:
                    # check the right child
                    if interval['end'] > self.nodelist[curr_node].subtreeMax:
                        # check if the current subtree max is greatr than our
                        # end, if not set it to our end
                        self.nodelist[curr_node].subtreeMax = interval['end']
                    if self.nodelist[curr_node].right == -1:
                        # if the right child link is empty, add this node
                        # get indicies, set right child link to equal
                        # this index
                        my_index = len(self.nodelist)
                        self.nodelist[curr_node].right = my_index
                        # add to tree
                        self.nodelist.append(
                                MyNode(interval, my_index, curr_min, curr_max)
                                )
                        # done! break out of loop
                        break
                    else:
                        # if right child was not empty, check the right subtree
                        curr_node = self.nodelist[curr_node].right

    def rebuild(self, sort_nodes, start_index, end_index):
        """recursive function used in balancing the tree. Input\
                sorted array of nodes, start and end index of that array
                to rebuild. Relinks nodes to produce an optimally balanced tree"""
        # check base case
        if start_index > end_index:
            return []

        # find mid index of our sub-array, set the current root to that index
        mid_index = (start_index + end_index) // 2
        curr_root = sort_nodes[mid_index]

        # get left values of the current root by calling the function
        # on the lower half of the sorted array
        left_values = self.rebuild(sort_nodes, start_index, mid_index - 1)
        # check if we returned a subtree of nodes and link it to our current root
        if left_values != []:
            # set out current root's left link to recursive call's root
            self.nodelist[curr_root].left = left_values[0] 

            # propogate subtreeMin value up tree
            # if our current min is smaller, propogate that instead
            if left_values[1] < self.nodelist[curr_root].start:
                self.nodelist[curr_root].subtreeMin = left_values[1] 
            else:
                self.nodelist[curr_root].subtreeMin = self.nodelist[curr_root].start 

        else:
            # if we didn't create a left subtree, set our left child link to -1
            # and propogate our subtreeMin
            self.nodelist[curr_root].left = -1 
            self.nodelist[curr_root].subtreeMin = self.nodelist[curr_root].start

        # get right child subtree with recursive call to upper half of array
        right_values = self.rebuild(sort_nodes, mid_index + 1, end_index)
        if right_values != []:
            # if we got a subtree, set our current root's right child to 
            # that subtree
            self.nodelist[curr_root].right = right_values[0]

            # propogate up subtreeMax if it is greater than our end, if not
            # propogate our end as the subtreeMax
            if right_values[2] > self.nodelist[curr_root].end:
                self.nodelist[curr_root].subtreeMax = right_values[2]
            else:
                self.nodelist[curr_root].subtreeMax = self.nodelist[curr_root].end
        else:
            # if we didn't get a right subtree, set our right child link to
            # -1 and propogate our end as subtreeMax
            self.nodelist[curr_root].right = -1
            self.nodelist[curr_root].subtreeMax = self.nodelist[curr_root].end

        # collect our current subtree min and max after all recursive calls
        min_val = self.nodelist[curr_root].subtreeMin
        max_val = self.nodelist[curr_root].subtreeMax
        # return this to previous function call along with our current root index
        return [curr_root, min_val, max_val]

    def inOrderTraverse(self):
        """ perform iterativein-order traversal of tree and return a sorted array\
            of tree node indicies """

        # initialize sort nodes list to empty
        sort_nodes = []
        if self.nodelist == []:
            # if no nodes in tree, return empty list
            return [] 
        else:
            # set current node to root, initialize the stack to empty
            curr_node = self.root
            node_stack = []
            done = False
            # iterate until we finish the traversal
            while not done:
                if curr_node != -1:
                    # if the current node exists, add to stack and 
                    # set curr_node to it's left child
                    node_stack.append(curr_node)
                    curr_node = self.nodelist[curr_node].left
                else:
                    # output, then check right node

                    # check length of stack: if > 0, then process last item
                    if len(node_stack) > 0:
                        # pop a node off the stack and append to the sorted
                        # node array
                        curr_node = node_stack.pop()
                        sort_nodes.append(curr_node)

                        # set current node to its right child
                        curr_node = self.nodelist[curr_node].right
                    else:
                        # if stack is empty, make this the last loop
                        done = True

        # return sorted array of nodes
        return sort_nodes    

    def balance(self):
        """function to balance the tree"""
        # perform in-order traversal and get a sorted list of nodes
        sort_nodes = self.inOrderTraverse()
        # if the tree is empty, return out of function
        if sort_nodes == []:
            return
        # now rebuild tree using recursive rebuild function
        root_values = self.rebuild(sort_nodes, 0,len(sort_nodes)-1)
        # set our root to root from initial recursive call
        self.root = root_values[0]
            
    def nodeOverlap(self, interval, curr_node):
        """function to check overlap between an interval and a node and\
                return a boolean"""
        # if there is overlap, return true, else return false
        if (self.nodelist[curr_node].start < interval['end'] and
                self.nodelist[curr_node].end > interval['start']):
            return True
        return False

## END BASE CLASSES ##

## FINAL CLASSES: Inherit from Base Classes ##
class ExonTree(MyTree):
    """ExonTree Object holds all exons for a gene"""
    def __init__(self):
        # Initialize 
        MyTree.__init__(self)

    def findNodeBool(self, interval):
        """Find if an interval overlaps with any nodes in this tree\
                , return a boolean whether it was found"""

        # if tree is empty, return False
        if self.nodelist == []:
            return False
        else:
            # set current node to root, initialize node stack
            curr_node = self.root
            node_stack = []
            done = False

            while not done:
                # while stack is nonempty, iterate
                if curr_node != -1:
                    # if current node is valid, add to stack
                    node_stack.append(curr_node)

                    # get left child index
                    left = self.nodelist[curr_node].left

                    # check if left subtree's max is greater than the interval start
                    if self.nodelist[left].subtreeMax > interval['start']:
                        # set current node to left child for next loop
                        curr_node = self.nodelist[curr_node].left
                    else:
                        # if not, set current node to invalid index
                        curr_node = -1
                else:
                    # if current node is an invalid index, check stack
                    if len(node_stack) > 0:
                        # if stack is nonempty, process last item on stack
                        curr_node = node_stack.pop()

                        # if this node DOES overlap with the interval,
                        # short-circuit find operation and return True
                        if self.nodeOverlap(interval, curr_node):
                            return True

                        # if not, grab right index
                        right = self.nodelist[curr_node].right
                        if (right != -1 and
                                self.nodelist[right].subtreeMin < interval['end']):
                            # if right child is valid and has a subtreeMIn that
                            # is less than current node's end,
                            # set current node to right child
                            curr_node = self.nodelist[curr_node].right
                        else:
                            # set current node to invalid index
                            curr_node = -1
                    else:
                        # if stack is empty, end loop
                        done = True

            # if we didn't short circuit, that means no nodes overlapped with
            # the interval, return False
            return False

class GeneNode(MyNode):
    """GeneNode class inherits from MyNode class and stores extra information\
            about the gene, as well as handles exons of the gene and implements\
            an addReads operation and output method"""

    def __init__(self,interval,geneid,gene_index,submin,submax,num_of_bams):
        """initialization function"""
        # initialize base class
        MyNode.__init__(self,interval,gene_index,submin,submax)
        # set gene data 
        self.gene_id = geneid
        # create an exon tree for this gene
        self.exons = ExonTree()
        # initialize read counts to 0 for each bam
        self.reads = [0] * num_of_bams

    def addReads(self, bam_num):
        """add reads to a specific bam file index"""
        self.reads[bam_num] += 1

    def checkExons(self, interval):
        """check whether the interval overlaps with any exon of this gene"""
        return self.exons.findNodeBool(interval)

    def writeOut(self, out_fp, mychrom, mystrand):
        """ output function goes through all bam files and outputs to file """
        # construct output string
        outstr = mychrom + "\t" + str(self.start) + "\t" + str(self.end) +\
                "\t" + mystrand + "\t" + self.gene_id 
        
        # add read counts
        for i in self.reads:
            outstr += "\t" + str(i)

        outstr += "\n"

        # write out
        out_fp.write(outstr)

class GeneTree(MyTree):
    """GeneTree object inherits from MyTree and adds a new find operation, \
            more information about genes, adds GeneNode objects and an output\
            method"""
    def __init__(self, chrom, strand):
        """initialization function"""
        # initialize base class
        MyTree.__init__(self)
        # set shared chrom/strand information for all genes in this tree
        self.chrom = chrom
        self.strand = strand

    def addExon(self, fields, gene_index):
        """add an exon to this gene"""
        self.nodelist[gene_index].exons.addNode(fields)

    def overlapInterval(self,interval, bam_num):
        """Overlap an interval with the gene nodes in this tree, checking\
                whether it overlaps multiple intervals and also checking\
                if it overlaps exons. Add reads to appropriate bam file index\ 
                if this is true"""

        # find the list of indicies of nodes that overlap with this interval
        indicies = self.findNode(interval)

        # if we found no overlaps, return out of this function
        if indicies == []: 
            return False

        # initialize exon overlap checks
        count_exon_overlaps = 0
        true_gene = -1

        # if we had multiple gene overlaps:
        if len(indicies) > 1:
            # check whether the interval overlaps exons for all those genes
            for i in indicies:
                if self.nodelist[i].checkExons(interval):
                    # increment exon overlaps and set
                    # true_gene to this gene
                    count_exon_overlaps += 1
                    true_gene = i

            # if the interval overlapped multiple genes' exons or none,
            # return False
            if count_exon_overlaps > 1 or count_exon_overlaps == 0:
                return False
            elif count_exon_overlaps == 1:
                # else if only one gene's exons overlapped the interval
                # increment the read count for this bam file
                self.nodelist[true_gene].addReads(bam_num)
                return True

        
        # if only one gene overlapped:
        if self.nodelist[indicies[0]].checkExons(interval): 
            # if exons of that one gene overlapped, then increment reads
            self.nodelist[indicies[0]].addReads(bam_num)
            return True

        return False

    def writeTree(self, out_fp):
        """output function to write out the entire tree"""
        # sort the nodes for output by start/end position
        # with an in-order traversal
        sort_nodes = self.inOrderTraverse()

        # for each node, output its information to file
        for i in sort_nodes: 
            self.nodelist[i].writeOut(out_fp, self.chrom, self.strand)

    def findNode(self, interval):
        """function to find GeneNodes that overlap with the interval of interest\
                and return a list of overlapping node indicies"""

        # if nodelist is empty, return empty list
        if self.nodelist == []:
            return []
        else:
            # initialize current node to root and stack and overlap list to empty
            curr_node = self.root
            overlap_list = []
            node_stack = []
            done = False

            # while not done iterating:
            while not done:
                if curr_node != -1:
                    # if current node is valid, add it to the stack
                    node_stack.append(curr_node)
                    # get left child index
                    left = self.nodelist[curr_node].left
                    # check if left child node is valid and if subtreeMax is greater
                    # than interval start
                    if (left != -1 and
                            self.nodelist[left].subtreeMax > interval['start']):
                        # set current node to left child for next loop
                        curr_node = self.nodelist[curr_node].left
                    else:
                        # set current node to invalid
                        curr_node = -1
                else:
                    # check length of stack
                    if len(node_stack) > 0:
                        # if stack is nonempty, process last item on stack
                        curr_node = node_stack.pop()
                        
                        # if this node overlaps with interval:
                        if self.nodeOverlap(interval, curr_node):
                            # add it to the overlap list
                            overlap_list.append(self.nodelist[curr_node].index)

                        # get right index
                        right = self.nodelist[curr_node].right

                        # if right child is valid and right subtree min is less
                        # than our end:
                        if (right != -1 and
                                self.nodelist[right].subtreeMin < interval['end']):
                            # set current node to right child
                            curr_node = self.nodelist[curr_node].right
                        else:
                            # set current node to invalid index
                            curr_node = -1
                    else:
                        # if stack is empty, end loop
                        done = True

            # return list of overlapping node indicies
            return overlap_list


    def addNode(self, interval, gene_id,bam_num):
        """Add a GeneNode to the GeneTree. Exact same function as ExonTree's\
                AddNode, except we add a GeneNode with extra information instead\
                Comments are the same as ExonTree's addNode so I will not repeat them"""
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
        """Balance this tree and all exon subtrees"""
        # balance this GeneTree
        self.balance()

        for gene in self.nodelist:
            # for each gene, balance all exon subtrees
            gene.exons.balance()
