

class Path:

    def __init__(self, start_node, end_node, n):
        self.start_node = start_node
        self.end_node = end_node
        self.n = n
        self.path_function = None


    def set_path_function(self, g, edge_list):
        if edge_list[0][0] != self.start_node or edge_list[-1][1] != self.end_node:
            print("Not a valid path for this connection. Does not start or finish at correct node")
            raise ValueError
        self.path_function = {}
        edges = g.get_edges()
        for edge in edges:
            self.path_function[int(edge[0]), int(edge[1])] = 0
        for edge in edge_list:
            if edge[0] > edge[1]:
                self.path_function[int(edge[1]),int(edge[0])] = 1
            else:
                self.path_function[int(edge[0]), int(edge[1])] = 1

    def is_cyclic(self):
        if self.path_function == None:
            print("Please set the value of the path before checking")
            raise ValueError
        else:
            traversed_nodes = []
            for edge in self.path_function:
                if edge[1] in traversed_nodes:
                    return True
                # if empty then need to add start node too
                if not traversed_nodes:
                    traversed_nodes.add(edge[0])
                traversed_nodes.add(edge[1])
            return False

    def is_too_long(self, g, length_property, l_switch, l_inf):
        loss = 0.0
        edges = g.get_edges([length_property])
        for source, target, length in edges:
            if source > target:
                loss += self.path_function[int(target),int(source)] * 0.2 * length
            else:
                loss += self.path_function[int(source), int(target)] * 0.2 * length
        loss += (sum(self.path_function.values()) - 1) * l_switch
        if loss >= l_inf:
            return True
        else:
            return False


    def get_distance(self, g, length_property, l_switch):
        loss = 0.0
        edges = g.get_edges([length_property])
        for source, target, length in edges:
            if source > target:
                loss += self.path_function[int(target), int(source)] * 0.2 * length
            else:
                loss += self.path_function[int(source), int(target)] * 0.2 * length
        loss += (sum(self.path_function.values()) - 1) * l_switch
        return loss / 0.2