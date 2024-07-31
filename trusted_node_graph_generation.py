from graph_tool.all import *
from random_graph_generation import SpecifiedTopologyGraph
from path import Path
import csv
import os
import numpy as np
import pandas as pd
import networkx as nx
from math import factorial
from itertools import islice
import time


# For generating the graph information needed to investigate the optimal topologies for trusted node networks.

def generate_fully_connected_graph(n, db_switch, box_size):
    graph = SpecifiedTopologyGraph()
    # will be fully connected if r = sqrt(2) box_size i.e. the entire box is covered by the circle (including the diagonal)
    # this occurs when (av_conn = 2 pi box_size) as r = sqrt(av_conn/ (pi * n)) * box_size
    graph.generate_random_mesh_graph(n, no_of_bobs = 0, no_of_bob_locations=0, dbswitch = db_switch, box_size = box_size, no_of_conns_av = n * 2 * np.pi)
    return graph

def get_capacity_edge(graph, dictionary, switched = False, db_switch = 1):
    edges_array = graph.get_edges(eprops=[graph.lengths_of_connections, graph.lengths_with_switch])
    edges_array = np.array(edges_array)
    capacities = {}
    for edge in edges_array:
        distance = edge[2]
        if switched == True:
            distance = distance + 5 * db_switch
        length = round(distance, 2)
        if length > 999:
            capacity = 0.0
        else:
            # from the look-up table
            capacity = dictionary["L" + str(length)]
        if capacity > 0.00000001:
            if (edge[0],edge[1]) in capacities.keys():
                capacities[edge[0], edge[1]].append((capacity, edge[2]))
            else:
                capacities[edge[0], edge[1]] = [(capacity, edge[2])]
    return capacities


def store_capacities(capacities, store_location, graph_id):
    for key in capacities.keys():
        dictionary_fieldnames = ["ID", "source_node", "target_node", "capacity", "distance"]
        values = {"ID": graph_id, "source_node": key[0], "target_node": key[1], "capacity" : capacities[key][0][0], "distance": capacities[key][0][1]}
        dictionary = [values]
        if os.path.isfile(store_location + '.csv'):
            with open(store_location + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writerows(dictionary)
        else:
            with open(store_location + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
                writer.writerows(dictionary)


def store_key_dict_full_set(graph, store_location, graph_id):
    node_array = graph.get_vertices()
    node_array = np.array(node_array)
    dictionary_fieldnames = ["ID", "source", "target", "key_value"]
    for node_1 in node_array:
        for node_2 in node_array:
            if node_1 < node_2:
                key_values = {"ID": graph_id, "source": node_1, "target": node_2, "key_value": 1}
                dictionary = [key_values]
                if os.path.isfile(store_location + '.csv'):
                    with open(store_location + '.csv', mode='a') as csv_file:
                        writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                        writer.writerows(dictionary)
                else:
                    with open(store_location + '.csv', mode='a') as csv_file:
                        writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                        writer.writeheader()
                        writer.writerows(dictionary)

def store_key_dict_random_set_user_nodes(graph, number_users, store_location, graph_id):
    node_array = graph.get_vertices()
    user_nodes = np.random.choice(node_array, size=number_users, replace=False)
    user_nodes = np.array(user_nodes)
    dictionary_fieldnames = ["ID", "source", "target", "key_value"]
    for node_1 in user_nodes:
        for node_2 in user_nodes:
            if node_1 < node_2:
                key_values = {"ID": graph_id, "source": node_1, "target": node_2, "key_value": 1}
                dictionary = [key_values]
                if os.path.isfile(store_location + '.csv'):
                    with open(store_location + '.csv', mode='a') as csv_file:
                        writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                        writer.writerows(dictionary)
                else:
                    with open(store_location + '.csv', mode='a') as csv_file:
                        writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                        writer.writeheader()
                        writer.writerows(dictionary)

def store_position_data(graph, vertex_store_location, edge_store_location, graph_id):
    dictionary_fieldnames = ["ID", "node", "xcoord", "ycoord"]
    for node in range(len(graph.vertex_properties["x_coord"].a)):
        values = {"ID": graph_id, "node": node, "xcoord": graph.vertex_properties["x_coord"].a[node], "ycoord": graph.vertex_properties["y_coord"].a[node]}
        dictionary = [values]
        if os.path.isfile(vertex_store_location + '.csv'):
            with open(vertex_store_location + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writerows(dictionary)
        else:
            with open(vertex_store_location + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
                writer.writerows(dictionary)
    dictionary_fieldnames = ["ID", "source_node", "target_node", "length"]#
    edges_array = graph.get_edges(eprops=[graph.lengths_of_connections])
    edges_array = np.array(edges_array)
    for edge in edges_array:
        values= {"ID": graph_id, "source_node" : edge[0], "target_node": edge[1], "length": edge[2]}
        dictionary = [values]
        if os.path.isfile(edge_store_location + '.csv'):
            with open(edge_store_location + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writerows(dictionary)
        else:
            with open(edge_store_location + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
                writer.writerows(dictionary)

if __name__ == "__main__":
    dictionary_bb84 = {}
    with open('rates_hotbob_bb84_20_eff.csv', mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            dictionary_bb84["L" + str(round(float(row["L"]), 2))] = float(row['rate'])
            line_count += 1
        print(f'Processed {line_count} lines.')
    j = 0
    for n in np.arange(8,10,2):
        for i in range(30):
            graph = generate_fully_connected_graph(n = n, db_switch= 0, box_size = 100)
            # parameter sweep here::::
            capacities = get_capacity_edge(graph.graph, dictionary_bb84, switched = False)
            store_capacities(capacities, store_location="5_trusted_nodes_capacities", graph_id = i+j*30)
            store_position_data(graph.graph, vertex_store_location = "5_trusted_nodes_positions", edge_store_location= "5_trusted_nodes_edge_data", graph_id = i+j*30)
            store_key_dict_full_set(graph.graph, store_location = "5_trusted_nodes_key_dict", graph_id = i+j*30)
            # for db_switch in np.arange(0.25,3,0.25):
            #     capacities_switched = get_capacity_edge(graph.graph, dictionary_bb84, switched=True, db_switch = db_switch)
            #     store_capacities(capacities_switched, store_location=f"4_trusted_nodes_capacities_switched_db_{int(db_switch * 100)}", graph_id=i + j * 5)
        j+= 1