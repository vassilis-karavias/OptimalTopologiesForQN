# OptimalTopologiesForQN
A repository for determining the optimal topology of TN QNs for a fixed given network cost.  
## Requirements
numpy: 1.20.1+  
graph-tool: 2.37+  
pandas: 1.5.2+  
scipy: 1.7.3+  
cplex: V20.1  
matplotlib: 3.6.2+   
networkx: 2.8.8+  

## Preprocessing

To generate the necessary csv files for the optimisation use the python file *trusted_node_graph_generation.py*. First you need to import the key rates from the key rate csv file into a dictionary of the form *{"L+distance_rounded_to_2_dp": rate}*. Some example key rate files are provided. To do this you can use the following code:  
*dictionary_bb84 = {}  
    with open('rates_hotbob_bb84_20_eff.csv', mode='r') as csv_file:  
        csv_reader = csv.DictReader(csv_file)  
        line_count = 0  
        for row in csv_reader:  
            if line_count == 0:   
                line_count += 1  
            dictionary_bb84["L" + str(round(float(row["L"]), 2))] = float(row['rate'])   
            line_count += 1*   
Then to generate random graphs use the *graph = generate_fully_connected_graph(n, db_switch, box_size)* method in the file. This will create a random fully connected graph of n nodes, using a switch loss of *db_switch* contained within a box of size *box_size* in km. Note that if no switches are needed then you can set *db_switch* to 0.  
To get the capacities of the edges in the graph, you can use the method:  
*capacities = get_capacity_edge(graph.graph, dictionary_bb84, switched)*  
where *graph.graph* is the graph-tool graph version of the topology graph *graph* obtained from the *generate_fully_connected_graph* method. *switched* is a boolean denoting whether the current graph uses switches or not.  
To store the necessary csv files, use:  
*store_capacities(capacities, store_location, graph_id)*  
*store_position_data(graph.graph, vertex_store_location, edge_store_location, graph_id)*  
*store_key_dict_full_set(graph.graph, store_location, graph_id)*  
The first method will store the capacities in the csv file named *store_location* with the columns *[ID, source_node, target_node, capacity, distance]*. *ID* is the *graph_id* used as a unique identifier of a graph - this allows multiple graphs to be stores in a single .csv file. *source_node, target_node* denotes the source and target of the current edge of the graph and *capacity* is the edge capacity for a single pair of devices. Finally *distance* is the distance between the edges in km. The second method will store the vertex data in a csv file named *vertex_store_location* and the edge data in a csv file named *edge_store_location*. The vertex file will have columns *[ID, node, xcoord, ycoord]* which labels the graph, node of the graph, xcoordinate of the node and ycoordinate of the node respectively. The edge file will have columns *[ID, source_node, target_node, length]*. Finally, the third method will store the key dictionary data (i.e. which nodes need keys generated between each other) in a csv file nmaed *store_location*. This csv file will be of the form *[ID, source, target, key_value]* where *source* is the source of the needed commodity, *target* is the sink of the needed commodity and *key_value* will be a relative value of the key rate needed - currently this is only implemented for *key_value = 1*. Note that all the store location names do not need the *.csv* ending.  

## Optimisation  

To perform the optimisation, the files *trusted_node_network.py* and *switched_fully_trusted_network.py* are provided. The former has the methods for the trusted node no-switching optimisation and the latter is the switching optimisation. To use these, first import the data from the .csv files generated in preprocessing using  
*key_dict, g = utils.import_problem_trusted_nodes_given_key_dict(node_file_path,  edge_file_path, key_dict_file_path)*  
The *node_file_path* is the *vertex_store_location* csv file previously generated. The *edge_file_path* is the capacity edge file with columns *[ID, source_node, target_node, capacity, distance]* and the *key_dict_file_path* is the key_dict csv file generated in the preprocessing. This will output a dictionary of *key_dicts* which are labelled based on the *ID*, i.e. *key_dicts = {ID: key_dict_current}* and a dictionary of networkx graphs *g = {ID: graph}*. You can then set up the non-switched problem by:  
*prob = cplex.Cplex()  
optimisation = trusted_node_network.Trusted_Node_Optimisation_No_Switching(prob, g[key], key_dict[key])  
sol_dict, prob, time_taken = optimisation.trusted_node_optimisation_program(time_limit,
                                                                                    cdevices,
                                                                                    c_dist,
                                                                                    c_conn_dist,
                                                                                    cnetwork,
                                                                                    Lambda, Nmax, Tij)*  
The first step is to set up a *Cplex* class. Then set up the *Trusted_Node_Optimisation_No_Switching* class. This takes as input the Cplex class, the graph for the current *ID* and the key_dict for the current *ID*. Then the method *trusted_node_optimisation_program* in the *Trusted_Node_Optimisation_No_Switching* class can be used to optimise and get the solution. *time_limit* is the maximum runtime of the optimisation, *cdevices* is the cost of the pair of devices, i.e. source and detector pair, *c_dist* is the cost of the fibre per km, *c_conn_dist* is the initial cost of laying out fibre per km, *cnetwork* is the maximum tolerable network cost, *Lambda* is the maximum number of fibres along an edge, *Nmax* is the maximum number of devices on an edge and *Tij* is the key requirement dict. *Tij* can either be a floating point value meaning all commodities in *key_dict* have equal weight of importance or *Tij* is a dictionary of values meaning the importance of a given commodity (i.e. how many more keys are needed for that commodity) scales as Tij[comm_1] / Tij[comm_2].  
For the switching optimisation problem, you first need to make the *key_dict* bidirectional using:  
*key_dict_new = switched_fully_trusted_network.make_key_dict_bidirectional(key_dict[key])*  
Then you can set up the problem the same as the no-switching version, but using the switching models:  
*prob = cplex.Cplex()  
optimisation = switched_fully_trusted_network.Switched_Trusted_Node_Optimisation_Constant_Switching_Ratio(prob, g[key], key_dict_new)  
sol_dict, prob = optimisation.trusted_node_optimisation_program(time_limit, csource, cdet,
                                                                    c_dist, c_conn_dist, Lambda, cgraph,
                                                                    f_switch, Nmax, Tij)*  
This time the source and detector costs are separate in  *csource* and *cdet*. *Nmax* now represents the maximum number of devices on a node and *cgraph* is the equivalent of *cnetwork*. Finally, *f_switch* is the fraction of time lost to calibration.  
There are several methods in *trusted_node_network.py* that can be used to investigate the various parameters.
                                    
