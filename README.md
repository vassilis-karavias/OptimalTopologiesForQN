# OptimalTopologiesForQN
A repository for determining the optimal topology of TN QNs for a fixed given network cost.
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
where *graph.graph* is the networkx graph version of the topology graph *graph* obtained from the *generate_fully_connected_graph* method. *switched* is a boolean denoting whether the current graph uses switches or not.  
To store the necessary csv files, use:  
*store_capacities(capacities, store_location, graph_id)*  
*store_position_data(graph.graph, vertex_store_location, edge_store_location, graph_id)*  
*store_key_dict_full_set(graph.graph, store_location, graph_id)*
The first method will store the capacities in the csv file named *store_location* with the columns *[ID, source_node, target_node, capacity, distance]*. *ID* is the *graph_id* used as a unique identifier of a graph - this allows multiple graphs to be stores in a single .csv file. *source_node, target_node* denotes the source and target of the current edge of the graph and *capacity* is the edge capacity for a single pair of devices. Finally *distance* is the distance between the edges in km. The second method will store the vertex data in a csv file named *vertex_store_location* and the edge data in a csv file named *edge_store_location*. The vertex file will have columns *[ID, node, xcoord, ycoord]* which labels the graph, node of the graph, xcoordinate of the node and ycoordinate of the node respectively. The edge file will have columns *[ID, source_node, target_node, length]*. Finally, the third method will store the key dictionary data (i.e. which nodes need keys generated between each other) in a csv file nmaed *store_location*. This csv file will be of the form *[ID, source, target, key_value]* where *source* is the source of the needed commodity, *target* is the sink of the needed commodity and *key_value* will be a relative value of the key rate needed - currently this is only implemented for *key_value = 1*. Note that all the store location names do not need the *.csv* ending.  

