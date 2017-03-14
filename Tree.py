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

    def addNode(self, interval, index):
        if self.nodelist == []:
            self.nodelist.append(MyNode(interval, index, interval['start'], interval['end']))
            self.root = 0
        else:
            curr_node = self.root
            curr_max = self.nodelist[self.root].subtreeMax
            curr_min = self.nodelist[self.root].subtreeMin

            while(1):

                if self.comp(curr_node, interval) == -1:
                    if interval['start'] < self.nodelist[curr_node].subtreeMin:
                        self.nodelist[curr_node].subtreeMin = interval['start']

                    if self.nodelist[curr_node].left == -1:
                        self.nodelist[curr_node].left = len(self.nodelist)
                        self.nodelist.append(MyNode(interval, index, curr_min, curr_max))
                        break
                    else:
                        curr_node = self.nodelist[curr_node].left
                else:
                    if interval['end'] > self.nodelist[curr_node].subtreeMax:
                        self.nodelist[curr_node].subtreeMax = interval['end']
                    if self.nodelist[curr_node].right == -1:

                        self.nodelist[curr_node].right = len(self.nodelist)
                        self.nodelist.append(MyNode(interval, index, curr_min, curr_max))
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
            self.nodelist[curr_root].subtreeMin = left_values[1] 
        else:
            self.nodelist[curr_root].left = -1 
            self.nodelist[curr_root].subtreeMin = self.nodelist[curr_root].start

        right_values = self.recursive_rebuild(sort_nodes, mid_index + 1, end_index)
        if right_values != []:
            self.nodelist[curr_root].right = right_values[0]
            self.nodelist[curr_root].subtreeMax = right_values[2]
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


    def nodeOverlap(self, interval, curr_node):
        if (self.nodelist[curr_node].start < interval['end'] and
                self.nodelist[curr_node].end > interval['start']):
            return True
        return False

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
                    if (self.nodelist[curr_node].left != -1 and
                            self.nodelist[self.nodelist[curr_node].left].subtreeMax > interval['start']):
                        curr_node = self.nodelist[curr_node].left
                    else:
                        curr_node = -1
                else:
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

                
