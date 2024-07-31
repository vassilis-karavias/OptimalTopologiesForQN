import os
import csv
import pandas as pd
import networkx as nx
import numpy as np

class Graph_Collector():

    def __init__(self):
        """
        Initialises class for Graph Collector, here will take a set of data of form [Graph_ID, node] and [Graph_ID,
        source_node, target_node] and collect each of the graphs into a single large graph by assigning different IDs
        to each different node
        """
        self.maps = {}


    def map_nodes_for_node_set(self, graph_id, node_graph):
        """
        returns the new node label for the mapped node based on input graph_id and node_id in the graph - if this is an
        unseen node then add the new node to the map dictionary that maps graph_id, node_id -> node_id of collective
        graph

        Parameters
        ----------
        graph_id : ID of the graph (int)
        node_graph : ID of node in the graph (int)

        Returns : ID of node in the new combined graph (int)
        -------

        """
        # look for entry in self.maps - if it exists return the value in maps
        if str(graph_id) + "," + str(int(node_graph)) in self.maps:
            return self.maps[str(graph_id) + "," + str(int(node_graph))]
        else:
            # if not in maps then add the entry with value the next biggest value to that in maps
            if not bool(self.maps):
                self.maps[str(graph_id) + "," + str(int(node_graph))] = 0
                return 0
            else:
                highest = max(self.maps.values())
                self.maps[str(graph_id)+ "," + str(int(node_graph))] = highest + 1
                return highest + 1

    def map_nodes_for_edges(self, graph_id, node_graph):
        """
        Returns the node ID in the new graph for each element in an edge
        Parameters
        ----------
        graph_id : ID of the graph (int)
        node_graph : ID of node in the graph (int)

        Returns : ID of node in the new combined graph (int)
        -------

        """
        # if in self.maps return the value
        if str(int(graph_id)) + "," + str(int(node_graph)) in self.maps:
            return self.maps[str(int(graph_id)) + "," + str(int(node_graph))]
        else:
            # if not in self.maps add new value to self.maps and return this
            print("Warning: Key not found - means there is a node in edges not in node set")
            return self.map_nodes_for_node_set(graph_id, node_graph)


def import_problem_trusted_nodes(node_file_path, edge_file_path):
    edge_data = pd.read_csv(edge_file_path)
    node_data = pd.read_csv(node_file_path)
    ## need to map keys to new values which are correctly labelled as the detector nodes should not be included
    possible_ids = edge_data["ID"].unique()
    graphs = {}
    key_dict = {}
    for id in possible_ids:
        node_data_current = node_data[node_data["ID"] == id].drop(["ID"], axis=1)
        capacity_values_set = edge_data[edge_data["ID"] == id].drop(["ID"], axis=1)
        capacity_values_set = capacity_values_set.astype({"source_node": int, "target_node": int})
        graph = nx.from_pandas_edgelist(capacity_values_set, "source_node", "target_node", ["capacity", "distance"])
        graph = graph.to_undirected()
        graph = graph.to_directed()
        node_attr = node_data_current.set_index("node").to_dict("index")
        nx.set_node_attributes(graph, node_attr)
        graphs[id] = graph
        key_dict_current = {}
        for i in graph.nodes:
            for j in graph.nodes:
                if i<j:
                    key_dict_current[(int(i),int(j))] = 1
        key_dict[id] = key_dict_current
    return key_dict, graphs


def import_problem_trusted_nodes_given_key_dict(node_file_path, edge_file_path, key_dict_file_path):
    edge_data = pd.read_csv(edge_file_path)
    node_data = pd.read_csv(node_file_path)
    key_dict = pd.read_csv(key_dict_file_path)
    ## need to map keys to new values which are correctly labelled as the detector nodes should not be included
    possible_ids = edge_data["ID"].unique()
    graphs = {}
    key_dict_dict = {}
    for id in possible_ids:
        node_data_current = node_data[node_data["ID"] == id].drop(["ID"], axis=1)
        capacity_values_set = edge_data[edge_data["ID"] == id].drop(["ID"], axis=1)
        capacity_values_set = capacity_values_set.astype({"source_node": int, "target_node": int})
        graph = nx.from_pandas_edgelist(capacity_values_set, "source_node", "target_node", ["capacity", "distance"])
        graph = graph.to_undirected()
        graph = graph.to_directed()
        node_attr = node_data_current.set_index("node").to_dict("index")
        nx.set_node_attributes(graph, node_attr)
        graphs[id] = graph
        key_dict_current = key_dict[key_dict["ID"] == id].drop(["ID"], axis =1)
        key_dict_current_dict = {}
        for index, row in key_dict_current.iterrows():
            key_dict_current_dict[(row["source"],row["target"])] = row["key_value"]
        key_dict_dict[id] = key_dict_current_dict
    return key_dict_dict, graphs


def import_node_edge_data(node_path_file):
    node_data = pd.read_csv(node_path_file)#

    possible_ids = node_data["ID"].unique()
    positions = {}
    for id in possible_ids:
        node_positions = node_data[node_data["ID"] == id].drop(["ID"], axis=1)
        for index, row in node_positions.iterrows():
            if id in positions.keys():
                positions[id][row[0]] = np.array([row[1],row[2]])
            else:
                positions[id] = {row[0]: np.array([row[1],row[2]])}
    return positions


