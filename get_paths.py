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


def generate_fully_connected_graph(n, db_switch, box_size):
    graph = SpecifiedTopologyGraph()
    # will be fully connected if r = sqrt(2) box_size i.e. the entire box is covered by the circle (including the diagonal)
    # this occurs when (av_conn = 2 pi box_size) as r = sqrt(av_conn/ (pi * n)) * box_size
    graph.generate_random_mesh_graph(n, no_of_bobs = 0, no_of_bob_locations=0, dbswitch = db_switch, box_size = box_size, no_of_conns_av = n * 2 * np.pi)
    return graph

def gttonx(graph, length_of_connections, length_with_switches):
    ## convert graphtool graph to networkx graph
    edges_array = graph.get_edges(eprops = [length_of_connections, length_with_switches])
    edges_array = np.array(edges_array)
    edge_list = pd.DataFrame(data = edges_array, columns = ["source", "target", "weight", "weight_with_switches"])
    graph = nx.from_pandas_edgelist(edge_list, "source", "target", edge_attr=["weight", "weight_with_switches"],
                                    create_using=nx.Graph())
    return graph

def get_longest_path(graph, db_switch, l_inf):
    edge_weight = graph.edges.data("weight", default=1)
    sorted_weight  = sorted(edge_weight, key=lambda item: item[2])
    k = 0
    loss = 0.0
    for source, target, weight in sorted_weight:
        loss += weight * 0.2 + db_switch
        if loss < l_inf:
            k += 1
            continue
        else:
            return k

def f(k, n):
    if k == 2:
        return 1/(factorial(n-3))
    else:
        return 1/factorial(n-(k+1)) + f(k-1,n)


def k_shortest_paths(G, source, target, k, weight=None):
    return islice(nx.shortest_simple_paths(G, source, target, weight=weight), k)


def convert_to_edge_list(path):
    j = None
    edge_path = []
    for i in path:
        if j == None:
            j = i
        else:
            edge_path.append([j,i])
            j = i
    return edge_path


def get_valid_paths(graph, db_switch, l_inf):
    graph_network_x = gttonx(graph.graph.g, graph.graph.lengths_of_connections, graph.graph.lengths_with_switch)
    k = get_longest_path(graph_network_x, db_switch, l_inf)
    print("k is " + str(k))
    paths = {}
    if k == 1:
        n_paths_max  = 1
    elif k == 0:
        print("no viable paths with non-zero capacity")
        raise ValueError
    else:
        if len(graph_network_x.nodes) <= 3:
            print("Graph needs more than 3 nodes")
            raise ValueError
        else:
            no_nodes = len(graph_network_x.nodes)
            n_paths_max = int(1+ factorial(no_nodes - 2) * f(k, no_nodes))
    print("n_paths_max is " + str(n_paths_max))
    nodes = sorted(graph_network_x.nodes)
    for i in range(len(nodes)):
        for j in range(i+1, len(nodes)):
            t_1 = time.time()
            vertex_1 = nodes[i]
            vertex_2 = nodes[j]
            n = 0
            paths_v1v2 = k_shortest_paths(graph_network_x, vertex_1, vertex_2, n_paths_max, weight = "weight_with_switches")
            print("time to find shortest paths: " + str(time.time() - t_1))
            for path in paths_v1v2:
                # if the length of a single path is larger than the k for the longest viable path that
                # allows l < l_inf than this path and all subsequent paths are definitely not valid

                # if len(path) > k:
                #     print("connection " + str(i) + "," + str(j) + "breaks because path is too long")
                #     break
                current_path = Path(start_node = vertex_1, end_node = vertex_2, n = n)
                path = convert_to_edge_list(path)
                current_path.set_path_function(graph.graph.g, path)
                is_too_long = current_path.is_too_long(graph.graph.g, length_property=graph.graph.lengths_of_connections,
                                                   l_switch=db_switch, l_inf=l_inf)
                is_valid = not is_too_long
                # if is too long then this path and all subsequent paths will be too long
                if is_valid:
                    if n == 0:
                        paths[vertex_1, vertex_2] = [current_path]
                    else:
                        paths[vertex_1, vertex_2].append(current_path)
                    n += 1
                else:
                    print("connection "+ str(i)  + "," + str(j) + "breaks because path has too high a loss")
                    break
            print("completed connection " + str(i)  + "," + str(j) + " in "  + str(time.time()-t_1))
    return paths


def get_valid_paths_test(graph, db_switch, l_inf):
    graph_network_x = gttonx(graph.graph.g, graph.graph.lengths_of_connections, graph.graph.lengths_with_switch)
    k = get_longest_path(graph_network_x, db_switch, l_inf)
    paths = {}
    if k == 1:
        n_paths_max = 1
    elif k == 0:
        print("no viable paths with non-zero capacity")
        raise ValueError
    else:
        if len(graph_network_x.nodes) <= 3:
            print("Graph needs more than 3 nodes")
            raise ValueError
        else:
            no_nodes = len(graph_network_x.nodes)
            n_paths_max = int(1 + factorial(no_nodes - 2) * f(k, no_nodes))
    nodes = sorted(graph_network_x.nodes)
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            t_1 = time.time()
            vertex_1 = nodes[i]
            vertex_2 = nodes[j]
            n = 0
            paths_v1v2 = k_shortest_paths(graph_network_x, vertex_1, vertex_2, n_paths_max,
                                          weight="weight_with_switches")
            for path in paths_v1v2:
                # if the length of a single path is larger than the k for the longest viable path that
                # allows l < l_inf than this path and all subsequent paths are definitely not valid


                current_path = Path(start_node=vertex_1, end_node=vertex_2, n=n)
                path = convert_to_edge_list(path)
                current_path.set_path_function(graph.graph.g, path)
                is_too_long = current_path.is_too_long(graph.graph.g,
                                                       length_property=graph.graph.lengths_of_connections,
                                                       l_switch=db_switch, l_inf=l_inf)
                is_valid = not is_too_long
                # if is too long then this path and all subsequent paths will be too long
                if is_valid:
                    if n == 0:
                        paths[vertex_1, vertex_2] = [current_path]
                    else:
                        paths[vertex_1, vertex_2].append(current_path)
                    n += 1
            print("completed connection without assumptions " + str(i) + "," + str(j) + " in " + str(time.time() - t_1))
    return paths

# def get_valid_paths(graph, db_switch, l_inf):
#     vertices = graph.graph.g.get_vertices()
#     paths = {}
#     for i in range(len(vertices)):
#         for j in range(i+1,len(vertices)):
#             vertex_1 = vertices[i]
#             vertex_2 = vertices[j]
#             n = 0
#             paths_v1v2 = all_paths(graph.graph.g, vertex_1, vertex_2, edges = True)
#             for path in paths_v1v2:
#                 current_path = Path(start_node = vertex_1, end_node = vertex_2, n = n)
#                 current_path.set_path_function(graph.graph.g, path)
#                 is_cyclic = current_path.is_cyclic()
#                 if is_cyclic:
#                     print("The path is cyclic")
#                 is_too_long = current_path.is_too_long(graph.graph.g, length_property = graph.graph.lengths_of_connections, l_switch = db_switch, l_inf = l_inf)
#                 is_valid = not is_cyclic or is_too_long
#                 if is_valid:
#                     if n == 0:
#                         paths[vertex_1, vertex_2] = [current_path]
#                     else:
#                         path[vertex_1, vertex_2].append(current_path)
#                     n += 1
#     return paths


def get_capacities_paths(graph, paths, db_switch, dictionary):
    capacities = {}
    for key in paths.keys():
        capacities[key] = []
        for path in paths[key]:
            length = path.get_distance(graph.graph.g, length_property = graph.graph.lengths_of_connections, l_switch = db_switch)
            length = round(length, 2)
            if length > 999:
                capacity = 0.0
            else:
                # from the look-up table
                capacity = dictionary["L" + str(length)]
            if capacity > 0.00000001:
                capacities[key].append(capacity)
    return capacities

def store_paths(paths, capacities, path_storage):
    for key in paths.keys():
        for i in range(len(paths[key])):
            path = paths[key][i]
            path_function = path.path_function
            dictionary_fieldnames = ["source_node", "target_node", "n", "capacity"] + [f"delta_{i}_{j}" for i,j in path_function.keys()]
            values = {"source_node": key[0], "target_node": key[1], "n" : path.n, "capacity": capacities[key][i]}
            for l,j in path_function.keys():
                values[f"delta_{l}_{j}"] = path_function[l,j]
            dictionary = [values]
            if os.path.isfile(path_storage + '.csv'):
                with open(path_storage + '.csv', mode='a') as csv_file:
                    writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                    writer.writerows(dictionary)
            else:
                with open(path_storage + '.csv', mode='a') as csv_file:
                    writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                    writer.writeheader()
                    writer.writerows(dictionary)


dictionary_bb84 = {}
with open('rates_hotbob_bb84.csv', mode='r') as csv_file:
    csv_reader = csv.DictReader(csv_file)
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            print(f'Column names are {", ".join(row)}')
            line_count += 1
        dictionary_bb84["L" + str(round(float(row["L"]), 2))] = float(row['rate'])
        line_count += 1
    print(f'Processed {line_count} lines.')


graph = generate_fully_connected_graph(n=20, db_switch = 1, box_size = 100)
paths = get_valid_paths(graph, db_switch = 1, l_inf = 26.6)
capacities = get_capacities_paths(graph, paths, db_switch= 1, dictionary= dictionary_bb84)
store_paths(paths, capacities, path_storage = "test_paths_20_nodes")
# paths = get_valid_paths_test(graph, db_switch = 1, l_inf = 26.6)
# capacities = get_capacities_paths(graph, paths, db_switch= 1, dictionary= dictionary_bb84)
# store_paths(paths, capacities, path_storage = "test_paths_no_assumptions")