import cplex
import networkx as nx
import utils
import numpy as np
import time
import csv
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import os
import switched_fully_trusted_network
import pandas as pd
### CONSTRAINTS
### IMPORT GRAPH
### Remember here key_dict[i,j]: i<j only



class Trusted_Node_Optimisation():

    def __init__(self, prob, g, key_dict):
        self.prob = prob
        self.g = g
        self.key_dict = key_dict
        super().__init__()


    def create_sol_dict(self):
        """
        Create a dictionary with the solution of the parameters
        """
        names = self.prob.variables.get_names()
        values = self.prob.solution.get_values()
        sol_dict={names[idx]: (values[idx]) for idx in range(self.prob.variables.get_num())}
        return sol_dict

    def split_sol_dict(self, sol_dict):
        """
        Split the solution dictionary into 2 dictionaries containing the flow variables only and the binary variables only
        Parameters
        ----------
        sol_dict : The solution dictionary containing solutions to the primary flow problem

        Returns
        -------

        """
        flow_dict = {}
        numbers_dict = {}
        lambda_dict = {}
        delta_dict = {}
        for key in sol_dict:
            # get all keys that are flow and add to dictionary
            if key[0] == "x" or key[0] == "z":
                flow_dict[key] = sol_dict[key]
            elif key[0] == "N":
                # get the keys that are binary 'on' 'off' and add to dictionary
                numbers_dict[key] = sol_dict[key]
            elif key[0] == "d":
                delta_dict[key] = sol_dict[key]
            else:
                # get the keys that represent the lambda_{i,j} - representing number of detectors and add to dictionary
                lambda_dict[key] = sol_dict[key]
        return flow_dict, numbers_dict, lambda_dict, delta_dict


    def log_optimal_solution_to_problem(self, save_file, graph_id):
        sol_dict = self.create_sol_dict()
        flow_dict, numbers_dict, lambda_dict = self.split_sol_dict(sol_dict)
        ## if file does not exist - we wish to store information of q_{i,j,d}^{m}, w_{i,j,d}^{m}, lambda_{d}^{m}, delta_{i,j,d}^{m}
        ## need to generate all possible values of i,j,d available. Need to think of the most appropriate way to store this data.
        dict = {"ID": graph_id}
        dict.update(flow_dict)
        dict.update(lambda_dict)
        dict.update(numbers_dict)
        dictionary = [dict]
        dictionary_fieldnames = list(dict.keys())
        with open(save_file + '.csv', mode='a') as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
            writer.writeheader()
            writer.writerows(dictionary)


    def log_optimal_solution_to_problem_column_format(self, save_file, graph_id):
        sol_dict = self.create_sol_dict()
        flow_dict, numbers_dict, lambda_dict = self.split_sol_dict(sol_dict)
        ## if file does not exist - we wish to store information of q_{i,j,d}^{m}, w_{i,j,d}^{m}, lambda_{d}^{m}, delta_{i,j,d}^{m}
        ## need to generate all possible values of i,j,d available. Need to think of the most appropriate way to store this data.
        dictionary = []
        dictionary_fieldnames = ["ID", "variable_name", "variable_value"]
        for key in flow_dict.keys():
            dictionary.append({"ID": graph_id, "variable_name": key, "variable_value": flow_dict[key]})
        for key in lambda_dict.keys():
            dictionary.append({"ID": graph_id, "variable_name": key, "variable_value": lambda_dict[key]})
        for key in numbers_dict.keys():
            dictionary.append({{"ID": graph_id, "variable_name": key, "variable_value": numbers_dict[key]}})
        with open(save_file + '.csv', mode='a') as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
            writer.writeheader()
            writer.writerows(dictionary)

    def conservation_of_flow_constraint(self):
        pass

    def no_flow_into_source_out_sink(self):
        pass

    def cost_graph_limit(self, *args, **kwargs):
        pass

    def limited_flow_across_edge(self):
        pass

    def add_fibre_on_constraint(self, Lambda):
        pass

    def fibre_limitation(self, no_channels):
        pass

    def maximise_lowest_capacity_connection_objective(self):
        pass


class Trusted_Node_Optimisation_No_Switching(Trusted_Node_Optimisation):

    def __init__(self, prob, g, key_dict):
        super().__init__(prob = prob, g = g, key_dict = key_dict)



    def conservation_of_flow_constraint(self):
        """
        add the conservation of flow constraint:
        \sum_{m \in N(n)} X_{(n,m)}^{k} - X_{(m,n)}^{k} = 0
        """
        variable_names = [f'x{i}_{j}_k{k[0]}_{k[1]}' for k in self.key_dict for i, j in list(self.g.edges)]
        self.prob.variables.add(names=variable_names, types=[self.prob.variables.type.continuous] * len(variable_names))
        for n in self.g.nodes:
            for k in self.key_dict:
                if k[1] != n and k[0] != n:
                    flow = []
                    coeffs = []
                    for m in self.g.neighbors(n):
                        flow.append(f"x{n}_{m}_k{k[0]}_{k[1]}")
                        coeffs.append(1)
                        flow.append(f"x{m}_{n}_k{k[0]}_{k[1]}")
                        coeffs.append(-1)
                    lin_expressions = [cplex.SparsePair(ind=flow, val=coeffs)]
                    self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["E"], rhs=[0.])



    def no_flow_into_source_out_sink(self):
        """
        add the constraint to prevent flow into source and out of sink:
        X_{(n,i)}^{k=(i,j)} = 0
        X_{(j,n)}^{k=(i,j)} = 0
        """
        for k in self.key_dict:
            for n in self.g.neighbors(k[1]):
                variables = [f"x{k[1]}_{n}_k{k[0]}_{k[1]}"]
                coeffs = [1]
                lin_expressions = [cplex.SparsePair(ind=variables, val=coeffs)]
                self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["E"], rhs=[0.])
            for n in self.g.neighbors(k[0]):
                variables = [f"x{n}_{k[0]}_k{k[0]}_{k[1]}"]
                coeffs = [1]
                lin_expressions = [cplex.SparsePair(ind=variables, val=coeffs)]
                self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["E"], rhs=[0.])


    def cost_graph_limit(self, cdevices, cfibre, cconn, cnetwork, Lambda, * args, **kwargs):
        """
        Add the constraint that limits the cost of the full graph:
        \sum_{(i,j) \in E} C_{i,j}^{devices}\lambda_{i,j} + C_{i,j}^{fibre}N_{i,j}^{fibre} \leq C_{network}
        """
        # cfibre is a dictionary with fibre costs based on the edge
        # NEED A FUNCTION TO GET CFIBRE
        lambda_variables = []
        n_fibres = []
        delta_variables = []
        for i,j in self.g.edges:
            if i < j:
                lambda_variables.append(f"lambda_{i}_{j}")
                n_fibres.append(f"N_{i}_{j}")
                delta_variables.append(f"delta_{i}_{j}")
        self.prob.variables.add(names=lambda_variables, types=[self.prob.variables.type.integer] * len(lambda_variables), ub = [Lambda] * len(lambda_variables))
        self.prob.variables.add(names=n_fibres, types=[self.prob.variables.type.integer] * len(n_fibres))
        self.prob.variables.add(names=delta_variables, types=[self.prob.variables.type.binary] * len(delta_variables))
        variables = []
        val = []
        for i,j in self.g.edges:
            if i < j:
                variables.extend([f"lambda_{i}_{j}", f"N_{i}_{j}", f"delta_{i}_{j}"])
                val.extend([cdevices, cfibre[i,j], cconn[i,j]])
        lin_expressions = [cplex.SparsePair(ind = variables, val = val)]
        self.prob.linear_constraints.add(lin_expr = lin_expressions, senses = ["L"], rhs = [float(cnetwork)])



    def limited_flow_across_edge(self):
        """
        Constrains the maximum flow along an edge dependent on the capacity if the devices across the edge
        \sum_{k \in K} X_{(i,j)}^{k} + X_{(j,i)}^{k} \leq \lambda_{i,j} c_{i,j}
        """
        for i,j in self.g.edges:
            if i < j:
                variables = []
                val = []
                for k in self.key_dict:
                    variables.extend([f"x{i}_{j}_k{k[0]}_{k[1]}", f"x{j}_{i}_k{k[0]}_{k[1]}"])
                    val.extend([1,1])
                variables.append(f"lambda_{i}_{j}")
                val.append(-int(self.g.edges[[i,j]]["capacity"]))
                lin_expressions = [cplex.SparsePair(ind = variables, val = val)]
                self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["L"], rhs=[0.])

    def fibre_limitation(self, no_channels = 1):
        """
        limits the number of devices across an edge to be at most the number of fibres across an edge:
        \lambda_{i,j} \leq N_{i,j}^{fibre}
        """
        for i,j in self.g.edges:
            if i < j:
                variables = [f"lambda_{i}_{j}", f"N_{i}_{j}"]
                val = [1/ no_channels, -1]
                lin_expressions = [cplex.SparsePair(ind=variables, val=val)]
                self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["L"], rhs=[0.])


    def add_fibre_on_constraint(self, Lambda):
        """
        adds the constraint that a fibre must be on to have any fibres N_{i,j}^{fibres} \leq Lambda \delta_{i,j}
        """
        for i,j in self.g.edges:
            if i < j:
                variables = [f"N_{i}_{j}", f"delta_{i}_{j}"]
                val = [1, -Lambda]
                lin_expressions = [cplex.SparsePair(ind=variables, val=val)]
                self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["L"], rhs=[0.])



    def maximise_throughput_objective(self):
        """
        Adds the maximisation of total throughput objective:
        maximise: \sum_{k \in K: i>k} \sum_{m \in N(j)} x^{k=(i,j)}_{(m,j)}

        """
        obj_vals = []
        for k in self.key_dict:
            for j in self.g.neighbors(k[1]):
                obj_vals.append((f"x{j}_{k[1]}_k{k[0]}_{k[1]}", 1))
        self.prob.objective.set_linear(obj_vals)
        self.prob.objective.set_sense(self.prob.objective.sense.maximize)


    def maximise_lowest_capacity_connection_objective(self, Tij):
        """
            Adds the maximise lowest capacity objective:
            maximise z
            subject to:
            z \leq \sum_{n \in N(j)} x_{(n,j)}^{k=(i,j)}
            """
        variable_names = ['z']
        self.prob.variables.add(names=variable_names, types=[self.prob.variables.type.continuous])
        for k in self.key_dict:
            variables = []
            val = []
            for n in self.g.neighbors(k[1]):
                if Tij == None:
                    variables.append(f"x{n}_{k[1]}_k{k[0]}_{k[1]}")
                    val.append(-1)
                elif Tij[k] > 0.000001:
                    variables.append(f"x{n}_{k[1]}_k{k[0]}_{k[1]}")
                    val.append(-1/Tij[k])
            variables.append("z")
            val.append(1)
            lin_expressions = [cplex.SparsePair(ind=variables, val=val)]
            self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["L"], rhs=[0.])
        obj_vals = [("z",1)]
        self.prob.objective.set_linear(obj_vals)
        self.prob.objective.set_sense(self.prob.objective.sense.maximize)

    def get_cfibre(self, c_dist, c_conn_dist):
        """
        get the cost of the fibre
        """
        cfibre = {}
        cconn = {}
        for i,j in self.g.edges:
            distance = self.g.edges[[i,j]]["distance"]
            cfibre[i,j] = distance * c_dist
            cconn[i,j] = distance * c_conn_dist
        return cfibre, cconn


    def trusted_node_optimisation_program(self, time_limit = 1e5, cdevices = 1.5, c_dist = 0.01, c_conn_dist = 0.1, cnetwork= 50, Lambda = 10, Nmax = 10, no_channels = 1, Tij = None):
        """
        set up and solve the problem for minimising the overall cost of the network
        """
        t_0 = time.time()
        print("Start Optimisation")
        self.conservation_of_flow_constraint()
        # add_capacity_constraint(prob, self.g, self.key_dict, Lambda=Lambda)
        self.no_flow_into_source_out_sink()
        cfibre, cconn = self.get_cfibre(c_dist, c_conn_dist=c_conn_dist)
        self.cost_graph_limit(cdevices, cfibre, cconn, cnetwork, Nmax)
        # add_minimise_trusted_nodes_objective(prob, self.g)
        self.limited_flow_across_edge()
        self.add_fibre_on_constraint(Lambda)
        self.fibre_limitation(no_channels)

        self.maximise_lowest_capacity_connection_objective(Tij)
        # maximise_throughput_objective(prob, self.g, self.key_dict)


        self.prob.write("test_multi.lp")
        self.prob.parameters.lpmethod.set(3)
        self.prob.parameters.mip.limits.cutpasses.set(1)
        self.prob.parameters.mip.strategy.probe.set(-1)
        self.prob.parameters.mip.strategy.variableselect.set(4)
        self.prob.parameters.mip.strategy.kappastats.set(1)
        self.prob.parameters.mip.tolerances.mipgap.set(float(0.01))
        # prob.parameters.simplex.limits.iterations = 50
        self.prob.parameters.timelimit.set(time_limit)
        print(self.prob.parameters.get_changed())
        t_1 = time.time()
        print("Time to set up problem: " + str(t_1 - t_0))
        self.prob.solve()
        t_2 = time.time()
        print("Time to solve problem: " + str(t_2 - t_1))

        print(f"The minimum Cost of Network: {self.prob.solution.get_objective_value()}")
        print(f"Number of Variables = {self.prob.variables.get_num()}")
        print(f"Number of Conditions = {self.prob.linear_constraints.get_num()}")
        sol_dict = self.create_sol_dict()
        flow_dict, binary_dict, lambda_dict, delta_dict = self.split_sol_dict(sol_dict)
        trusted_nodes = 0
        for key in delta_dict:
            trusted_nodes += delta_dict[key]
        print(f"Number of Trusted Nodes = {trusted_nodes}")
        return sol_dict, self.prob, t_2-t_1


### Graph Drawing #####


def plot_graph(graph, delta_dict, save_extra = ""):
    graph = graph.to_undirected()
    pos = {}
    for node in graph.nodes:
        pos[node] = [graph.nodes[node]["xcoord"], graph.nodes[node]["ycoord"]]
    plt.figure()

    nx.draw_networkx_nodes(graph, pos, node_shape="o")
    on_edges =[]
    off_edges =[]
    #### consider this but for cold and hot.....
    for key in delta_dict:
        current_edge = key[6:].split("_")
        current_edge = (int(current_edge[0]), int(current_edge[1]))
        on_off = delta_dict[key]
        if on_off:
            on_edges.append(current_edge)
        else:
            off_edges.append(current_edge)
    nx.draw_networkx_edges(graph, pos, edgelist=on_edges, edge_color="r", label = "On Edges")
    nx.draw_networkx_edges(graph, pos, edgelist=off_edges, edge_color="k", label = "Off Edges")
    plt.axis("off")
    plt.legend(loc = "best", fontsize = "small")
    plt.savefig(f"plot_graph_{save_extra}.jpg")
    plt.show()


### Investigation of parameters
### Parameters include: T_{i,j} and network costs

def cost_investigations(c_source, c_detector, c_fibre_per_km, c_initial_fibre_conn_per_km, overall_cost, Lambda, Tij, node_file_path, edge_file_path, key_dict_file_path, data_storage_location_keep_each_loop):
    key_dict, g = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path,
                                                                    edge_file_path=edge_file_path,
                                                                    key_dict_file_path=key_dict_file_path)
    if data_storage_location_keep_each_loop != None:
        if os.path.isfile(data_storage_location_keep_each_loop + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_detector_cost = last_row_explored["Cost_Detector"].iloc[0]
            current_fibre_per_km = last_row_explored["Cost_fibre_per_km"].iloc[0]
            current_init_fibre_cost = last_row_explored["C_Init_Fibre_conn_per_km"].iloc[0]
            current_key = last_row_explored["Graph key"].iloc[0]
        else:
            current_key = None
            current_detector_cost = 0.0
            current_fibre_per_km = 0.0
            current_init_fibre_cost = 0.0
            dictionary_fieldnames = ["Cost_Detector", "Cost_fibre_per_km", "C_Init_Fibre_conn_per_km", "Graph key", "objective_value"]
            with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_key = None
        current_detector_cost = 0.0
        current_fibre_per_km = 0.0
        current_init_fibre_cost = 0.0

    no_solution_list = []
    objective_list = {}
    for c_det in np.arange(current_detector_cost, c_detector, 0.1):
        for c_fibre in np.arange(current_fibre_per_km, c_fibre_per_km, 0.025):
            for c_initial_fibre in np.arange(current_init_fibre_cost, c_initial_fibre_conn_per_km, 0.05):
                for key in g.keys():
                    if current_key != None and current_key != key:
                        continue
                    elif current_key != None:
                        current_key = None
                        continue
                    if key not in no_solution_list:
                        try:
                            prob = cplex.Cplex()
                            optimisation = Trusted_Node_Optimisation_No_Switching(prob, g[key], key_dict[key])
                            sol_dict, prob, time_taken = optimisation.trusted_node_optimisation_program(time_limit = 1e3, cdevices = c_source + c_det, c_dist = c_fibre, c_conn_dist = c_initial_fibre, cnetwork= overall_cost, Lambda = Lambda, Nmax = Lambda, Tij = Tij)
                            flow_dict, binary_dict, lambda_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                            objective_value = prob.solution.get_objective_value()
                            ### step 2.
                            if key == 2:
                                plot_graph(g[key], delta_dict=delta_dict, save_extra=f"{c_det}_{c_fibre}_{c_initial_fibre}")
                            if c_det not in objective_list.keys():
                                objective_list[c_det] = {c_fibre: {c_initial_fibre: {key: objective_value}}}
                            elif c_fibre not in objective_list[c_det].keys():
                                objective_list[c_det][c_fibre] = {c_initial_fibre: {key: objective_value}}
                            elif c_initial_fibre not in objective_list[c_det][c_fibre].keys():
                                objective_list[c_det][c_fibre][c_initial_fibre] = {key: objective_value}
                            else:
                                objective_list[c_det][c_fibre][c_initial_fibre][key] = objective_value
                            if data_storage_location_keep_each_loop != None:
                                dictionary = [
                                    {"Cost_Detector" : c_det, "Cost_fibre_per_km": c_fibre, "C_Init_Fibre_conn_per_km": c_initial_fibre, "Graph key": key, "objective_value": objective_value}]
                                dictionary_fieldnames = ["Cost_Detector", "Cost_fibre_per_km", "C_Init_Fibre_conn_per_km", "Graph key", "objective_value"]
                                if os.path.isfile(data_storage_location_keep_each_loop + '.csv'):
                                    with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                                        writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                        writer.writerows(dictionary)
                                else:
                                    with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                                        writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                        writer.writeheader()
                                        writer.writerows(dictionary)
                        except:
                            no_solution_list.append(key)
                            continue
    if data_storage_location_keep_each_loop != None:
        plot_information = pd.read_csv(data_storage_location_keep_each_loop + ".csv")
        for index, row in plot_information.iterrows():
            if row["Cost_Detector"] not in objective_list.keys():
                objective_list[row["Cost_Detector"]] = {row["Cost_fibre_per_km"]: {row["C_Init_Fibre_conn_per_km"]: {row["Graph key"]: row["objective_value"]}}}
            elif row["Cost_fibre_per_km"] not in objective_list[row["Cost_Detector"]].keys():
                objective_list[row["Cost_Detector"]][row["Cost_fibre_per_km"]] = {row["C_Init_Fibre_conn_per_km"]: {row["Graph key"]: row["objective_value"]}}
            elif c_initial_fibre not in objective_list[row["Cost_Detector"]][row["Cost_fibre_per_km"]].keys():
                objective_list[row["Cost_Detector"]][row["Cost_fibre_per_km"]][row["C_Init_Fibre_conn_per_km"]] = {row["Graph key"]: row["objective_value"]}
            else:
                objective_list[row["Cost_Detector"]][row["Cost_fibre_per_km"]][row["C_Init_Fibre_conn_per_km"]][row["Graph key"]] = row["objective_value"]
        #### figure out what investigation to do here
        #### Firstly- this is relative to the cost of the source... c_{det}, c_{fibre_per_km}, c_{init_fibre}
        #### 1. Look at how each of these terms variation affects the overall graph cost - which one has the biggest effect
        #### 2. Investigate relationships between the cost terms, e.g. if the fibre placement cost scales the optimal
        #### solution will inevitably result in different graphs. How do the structures of these graphs vary with this cost
        #### And same for detector cost
        #### 3. Plot a set of 3-D graphs of max min capacity based on 2 parameters...

        #### Standard price for source: 1.5k, Standard price for detectors 120k, per fibre km: 3.8k, init cost to bury: 32.7k
        #### if we set the price of source to 0.01 --> det price: 0.8, per fibre km: 0.025, init cost to bury: 0.22 (set this to 0.2)

        #### Program for 1.
        objectives_detector_cost = {}
        objective_costs_to_normalise  = {}
        for cost in objective_list.keys():
            for key in objective_list[cost][0.025][0.2].keys():
                if cost < 0.00001:
                    objective_costs_to_normalise[key] = objective_list[cost][0.025][0.2][key]
                if objective_costs_to_normalise[key] > 0.00001:
                    if cost not in objectives_detector_cost.keys():
                        objectives_detector_cost[cost] = [objective_list[cost][0.025][0.2][key]/objective_costs_to_normalise[key]]
                    else:
                        objectives_detector_cost[cost].extend([objective_list[cost][0.025][0.2][key]/objective_costs_to_normalise[key]])
        mean_values = []
        std_values = []
        x = []
        for cost in objective_costs_to_normalise.keys():
            mean_values.append(np.mean(objectives_detector_cost[cost]))
            std_values.append(np.std(objectives_detector_cost[cost]))
            x.append(cost)

        plt.errorbar(x, mean_values, yerr=std_values, color="r", capsize=10)
        plt.xlabel("Cost of Detectors", fontsize=10)
        plt.ylabel("z_min (bits/s)", fontsize=10)
        plt.savefig("detector_cost_investigation.png")
        plt.show()

        objectives_fibre_cost = {}
        objective_costs_to_normalise = {}
        for cost in objective_list[0.8].keys():
            for key in objective_list[0.8][cost][0.2].keys():
                if cost < 0.00001:
                    objective_costs_to_normalise[key] = objective_list[0.8][cost][0.2][key]
                if objective_costs_to_normalise[key] > 0.00001:
                    if cost not in objectives_fibre_cost.keys():
                        objectives_fibre_cost[cost] = [
                            objective_list[0.8][cost][0.2][key] / objective_costs_to_normalise[key]]
                    else:
                        objectives_fibre_cost[cost].extend(
                            [objective_list[0.8][cost][0.2][key] / objective_costs_to_normalise[key]])
        mean_values = []
        std_values = []
        x = []
        for cost in objective_costs_to_normalise.keys():
            mean_values.append(np.mean(objectives_fibre_cost[cost]))
            std_values.append(np.std(objectives_fibre_cost[cost]))
            x.append(cost)

        plt.errorbar(x, mean_values, yerr=std_values, color="r", capsize=10)
        plt.xlabel("Cost of Fibre per km", fontsize=10)
        plt.ylabel("z_min (bits/s)", fontsize=10)
        plt.savefig("fibre_cost_investigation.png")
        plt.show()



        objectives_init_fibre_cost = {}
        objective_costs_to_normalise = {}
        for cost in objective_list[0.8][0.025].keys():
            for key in objective_list[0.8][0.025][cost].keys():
                if cost < 0.00001:
                    objective_costs_to_normalise[key] = objective_list[0.8][0.025][cost][key]
                if objective_costs_to_normalise[key] > 0.00001:
                    if cost not in objectives_init_fibre_cost.keys():
                        objectives_init_fibre_cost[cost] = [
                            objective_list[0.8][0.025][cost][key] / objective_costs_to_normalise[key]]
                    else:
                        objectives_init_fibre_cost[cost].extend(
                            [objective_list[0.8][0.025][cost][key] / objective_costs_to_normalise[key]])
        mean_values = []
        std_values = []
        x = []
        for cost in objective_costs_to_normalise.keys():
            mean_values.append(np.mean(objectives_init_fibre_cost[cost]))
            std_values.append(np.std(objectives_init_fibre_cost[cost]))
            x.append(cost)

        plt.errorbar(x, mean_values, yerr=std_values, color="r", capsize=10)
        plt.xlabel("Initial Cost of Fibre per km", fontsize=10)
        plt.ylabel("z_min (bits/s)", fontsize=10)
        plt.savefig("init_fibre_cost_investigation.png")
        plt.show()


        #### step 3.

        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        X = np.arange(0.0, c_fibre_per_km, 0.025)
        Y = np.arange(0.0, c_initial_fibre_conn_per_km, 0.05)
        X,Y = np.meshgrid(X,Y)
        Z = np.empty((len(X),len(Y)))
        for x in range(len(X)):
            for y in range(len(Y)):
                Z[x][y] = np.std(list(objective_list[0.8][X[x]][Y[y]].values()))
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        zmin, zmax = Z.min(), Z.max()
        sig_fig_min = int(np.floor(np.log10(abs(zmin))))
        sig_fig_max = int(np.floor(np.log10(abs(zmax))))
        ax.set_zlim(np.floor(zmin/sig_fig_min) * sig_fig_min, np.ceil(zmax /sig_fig_max) * sig_fig_max)
        ax.zaxis.set_major_locator(LinearLocator(10))
        # A StrMethodFormatter is used automatically
        ax.zaxis.set_major_formatter('{x:.02f}')
        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5)
        ax.set(
            xlabel='Fibre cost per km',
            ylabel='Initial fibre cost per km',
            zlabel='z_min (bits/s)',
        )
        plt.savefig("3d_plot_fibres.png")
        plt.show()


        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        X = np.arange(0.0, c_fibre_per_km, 0.025)
        Y = np.arange(0.0, c_detector, 0.1)
        X,Y = np.meshgrid(X,Y)
        Z = np.empty((len(X),len(Y)))
        for x in range(len(X)):
            for y in range(len(Y)):
                Z[x][y] = np.std(list(objective_list[Y[y]][X[x]][0.2].values()))
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        zmin, zmax = Z.min(), Z.max()
        sig_fig_min = int(np.floor(np.log10(abs(zmin))))
        sig_fig_max = int(np.floor(np.log10(abs(zmax))))
        ax.set_zlim(np.floor(zmin/sig_fig_min) * sig_fig_min, np.ceil(zmax /sig_fig_max) * sig_fig_max)
        ax.zaxis.set_major_locator(LinearLocator(10))
        # A StrMethodFormatter is used automatically
        ax.zaxis.set_major_formatter('{x:.02f}')
        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5)
        ax.set(
            xlabel='Fibre cost per km',
            ylabel='Detector cost',
            zlabel='z_min (bits/s)',
        )
        plt.savefig("3d_plot_fibres_detectors.png")
        plt.show()





        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        X = np.arange(0.0, c_initial_fibre_conn_per_km, 0.05)
        Y = np.arange(0.0, c_detector, 0.1)
        X,Y = np.meshgrid(X,Y)
        Z = np.empty((len(X),len(Y)))
        for x in range(len(X)):
            for y in range(len(Y)):
                Z[x][y] = np.std(list(objective_list[Y[y]][0.025][X[x]].values()))
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        zmin, zmax = Z.min(), Z.max()
        sig_fig_min = int(np.floor(np.log10(abs(zmin))))
        sig_fig_max = int(np.floor(np.log10(abs(zmax))))
        ax.set_zlim(np.floor(zmin/sig_fig_min) * sig_fig_min, np.ceil(zmax /sig_fig_max) * sig_fig_max)
        ax.zaxis.set_major_locator(LinearLocator(10))
        # A StrMethodFormatter is used automatically
        ax.zaxis.set_major_formatter('{x:.02f}')
        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5)
        ax.set(
            xlabel='Initial fibre cost per km',
            ylabel='Detector cost',
            zlabel='z_min (bits/s)',
        )
        plt.savefig("3d_plot_init_fibres_detectors.png")
        plt.show()


def cost_investigations_shorter(c_source, c_detector, c_fibre_per_km, c_initial_fibre_conn_per_km, overall_cost, Lambda, Tij, node_file_path, edge_file_path, key_dict_file_path, data_storage_location_keep_each_loop_fibre_init_fibre,data_storage_location_keep_each_loop_detect_init_fibre,data_storage_location_keep_each_loop_detect_fibre ):
    key_dict, g = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path,
                                                                    edge_file_path=edge_file_path,
                                                                    key_dict_file_path=key_dict_file_path)
    if data_storage_location_keep_each_loop_fibre_init_fibre != None:
        if os.path.isfile(data_storage_location_keep_each_loop_fibre_init_fibre + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop_fibre_init_fibre + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_fibre_per_km = last_row_explored["Cost_fibre_per_km"].iloc[0]
            current_init_fibre_cost = last_row_explored["C_Init_Fibre_conn_per_km"].iloc[0]
            current_key = last_row_explored["Graph key"].iloc[0]
        else:
            current_key = None
            current_fibre_per_km = 0.0
            current_init_fibre_cost = 0.0
            dictionary_fieldnames = ["Cost_fibre_per_km", "C_Init_Fibre_conn_per_km", "Graph key", "objective_value"]
            with open(data_storage_location_keep_each_loop_fibre_init_fibre + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_key = None
        current_fibre_per_km = 0.0
        current_init_fibre_cost = 0.0
    # for c_det in np.arange(current_detector_cost, c_detector, 0.1):
    no_solution_list = []
    objective_list = {}

    for c_fibre in np.arange(current_fibre_per_km, c_fibre_per_km, 0.025):
        for c_initial_fibre in np.arange(0.0, c_initial_fibre_conn_per_km, 0.05):
            if np.abs(c_fibre - current_fibre_per_km) < 0.0001:
                if c_initial_fibre < current_init_fibre_cost - 0.0001:
                    continue
            for key in g.keys():
                if current_key != None and current_key != key:
                    continue
                elif current_key != None:
                    current_key = None
                    continue
                if key not in no_solution_list:
                    try:
                        prob = cplex.Cplex()
                        optimisation = Trusted_Node_Optimisation_No_Switching(prob, g[key], key_dict[key])
                        sol_dict, prob, time_taken = optimisation.trusted_node_optimisation_program(time_limit = 2e2, cdevices = c_source + 0.8, c_dist = c_fibre, c_conn_dist = c_initial_fibre, cnetwork= overall_cost, Lambda = Lambda, Nmax = Lambda, Tij = Tij)
                        flow_dict, binary_dict, lambda_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                        objective_value = prob.solution.get_objective_value()
                        ### step 2.
                        # if key == 2:
                        #     plot_graph(g[key], delta_dict=delta_dict, save_extra=f"{c_fibre}_{c_initial_fibre}")
                        if c_fibre not in objective_list.keys():
                            objective_list[c_fibre] = {c_initial_fibre: {key: objective_value}}
                        elif c_initial_fibre not in objective_list[c_fibre].keys():
                            objective_list[c_fibre][c_initial_fibre] = {key: objective_value}
                        else:
                            objective_list[c_fibre][c_initial_fibre][key] = objective_value
                        if data_storage_location_keep_each_loop_fibre_init_fibre != None:
                            dictionary = [
                                {"Cost_fibre_per_km": c_fibre, "C_Init_Fibre_conn_per_km": c_initial_fibre, "Graph key": key, "objective_value": objective_value}]
                            dictionary_fieldnames = ["Cost_fibre_per_km", "C_Init_Fibre_conn_per_km", "Graph key", "objective_value"]
                            if os.path.isfile(data_storage_location_keep_each_loop_fibre_init_fibre + '.csv'):
                                with open(data_storage_location_keep_each_loop_fibre_init_fibre + '.csv', mode='a') as csv_file:
                                    writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                    writer.writerows(dictionary)
                            else:
                                with open(data_storage_location_keep_each_loop_fibre_init_fibre + '.csv', mode='a') as csv_file:
                                    writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                    writer.writeheader()
                                    writer.writerows(dictionary)
                    except:
                        no_solution_list.append(key)
                        continue
    if data_storage_location_keep_each_loop_fibre_init_fibre != None:
        plot_information = pd.read_csv(data_storage_location_keep_each_loop_fibre_init_fibre + ".csv")
        for index, row in plot_information.iterrows():
            if round(row["Cost_fibre_per_km"],3) not in objective_list.keys():
                objective_list[round(row["Cost_fibre_per_km"],3)] = {round(row["C_Init_Fibre_conn_per_km"],2): {row["Graph key"]: row["objective_value"]}}
            elif round(row["C_Init_Fibre_conn_per_km"],2) not in objective_list[round(row["Cost_fibre_per_km"],3)].keys():
                objective_list[round(row["Cost_fibre_per_km"],3)][round(row["C_Init_Fibre_conn_per_km"],2)] = {row["Graph key"]: row["objective_value"]}
            else:
                objective_list[round(row["Cost_fibre_per_km"],3)][round(row["C_Init_Fibre_conn_per_km"],2)][row["Graph key"]] = row["objective_value"]
        #### figure out what investigation to do here
        #### Firstly- this is relative to the cost of the source... c_{det}, c_{fibre_per_km}, c_{init_fibre}
        #### 1. Look at how each of these terms variation affects the overall graph cost - which one has the biggest effect
        #### 2. Investigate relationships between the cost terms, e.g. if the fibre placement cost scales the optimal
        #### solution will inevitably result in different graphs. How do the structures of these graphs vary with this cost
        #### And same for detector cost
        #### 3. Plot a set of 3-D graphs of max min capacity based on 2 parameters...

        #### Standard price for source: 1.5k, Standard price for detectors 120k, per fibre km: 3.8k, init cost to bury: 32.7k
        #### if we set the price of source to 0.01 --> det price: 0.8, per fibre km: 0.025, init cost to bury: 0.22 (set this to 0.2)

        #### Program for 1.
        # objectives_detector_cost = {}
        # objective_costs_to_normalise  = {}
        # for cost in objective_list.keys():
        #     for key in objective_list[cost][0.025][0.2].keys():
        #         if cost < 0.00001:
        #             objective_costs_to_normalise[key] = objective_list[cost][0.025][0.2][key]
        #         if objective_costs_to_normalise[key] > 0.00001:
        #             if cost not in objectives_detector_cost.keys():
        #                 objectives_detector_cost[cost] = [objective_list[cost][0.025][0.2][key]/objective_costs_to_normalise[key]]
        #             else:
        #                 objectives_detector_cost[cost].extend([objective_list[cost][0.025][0.2][key]/objective_costs_to_normalise[key]])
        # mean_values = []
        # std_values = []
        # x = []
        # for cost in objective_costs_to_normalise.keys():
        #     mean_values.append(np.mean(objectives_detector_cost[cost]))
        #     std_values.append(np.std(objectives_detector_cost[cost]))
        #     x.append(cost)
        #
        # plt.errorbar(x, mean_values, yerr=std_values, color="r", capsize=10)
        # plt.xlabel("Cost of Detectors", fontsize=10)
        # plt.ylabel("z_min (bits/s)", fontsize=10)
        # plt.savefig("detector_cost_investigation.png")
        # plt.show()

        objectives_fibre_cost = {}
        objective_costs_to_normalise = {}
        for cost in objective_list.keys():
            for key in objective_list[cost][0.2].keys():
                if cost < 0.00001:
                    objective_costs_to_normalise[key] = objective_list[cost][0.2][key]
                if objective_costs_to_normalise[key] > 0.00001:
                    if cost not in objectives_fibre_cost.keys():
                        objectives_fibre_cost[cost] = [
                            objective_list[cost][0.2][key] / objective_costs_to_normalise[key]]
                    else:
                        objectives_fibre_cost[cost].extend(
                            [objective_list[cost][0.2][key] / objective_costs_to_normalise[key]])
        mean_values = []
        std_values = []
        x = []
        for cost in objectives_fibre_cost.keys():
            mean_values.append(np.mean(objectives_fibre_cost[cost]))
            std_values.append(np.std(objectives_fibre_cost[cost]))
            x.append(cost)

        plt.errorbar(x, mean_values, yerr=std_values, color="r", capsize=10)
        plt.xlabel("Cost of Fibre per km", fontsize=10)
        plt.ylabel("z_min / z_min at zero cost", fontsize=10)
        plt.savefig("fibre_cost_investigation.png")
        plt.show()



        objectives_init_fibre_cost = {}
        objective_costs_to_normalise = {}
        for cost in objective_list[0.025].keys():
            for key in objective_list[0.025][cost].keys():
                if cost < 0.00001:
                    objective_costs_to_normalise[key] = objective_list[0.025][cost][key]
                if objective_costs_to_normalise[key] > 0.00001:
                    if cost not in objectives_init_fibre_cost.keys():
                        objectives_init_fibre_cost[cost] = [
                            objective_list[0.025][cost][key] / objective_costs_to_normalise[key]]
                    else:
                        objectives_init_fibre_cost[cost].extend(
                            [objective_list[0.025][cost][key] / objective_costs_to_normalise[key]])
        mean_values = []
        std_values = []
        x = []
        for cost in objectives_init_fibre_cost.keys():
            mean_values.append(np.mean(objectives_init_fibre_cost[cost]))
            std_values.append(np.std(objectives_init_fibre_cost[cost]))
            x.append(cost)

        plt.errorbar(x, mean_values, yerr=std_values, color="r", capsize=10)
        plt.xlabel("Initial Cost of Fibre per km", fontsize=10)
        plt.ylabel("z_min / z_min at zero cost", fontsize=10)
        plt.savefig("init_fibre_cost_investigation.png")
        plt.show()


        #### step 3.

        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        X = np.arange(0.0, c_fibre_per_km, 0.025)
        Y = np.arange(0.0, c_initial_fibre_conn_per_km, 0.05)
        X,Y = np.meshgrid(X,Y)
        Z = np.empty((len(X),len(X[0])))
        for i in range(len(X)):
            for j in range(len(X[i])):
                Z[i][j] = np.mean(list(objective_list[round(X[i][j],3)][round(Y[i][j],2)].values()))
        surf = ax.plot_surface(X, Y, Z, cmap = cm.coolwarm,
                               linewidth=0, antialiased=False)
        #
        zmin, zmax = Z.min(), Z.max()
        sig_fig_min = int(np.floor(np.log10(abs(zmin))))
        sig_fig_max = int(np.floor(np.log10(abs(zmax))))
        ax.set_zlim(np.floor(zmin/(10 ** sig_fig_min)) * (10 ** sig_fig_min), np.ceil(zmax/(10 ** (sig_fig_max-1))) * (10 ** (sig_fig_max-1)))
        ax.zaxis.set_major_locator(LinearLocator(10))
        # A StrMethodFormatter is used automatically
        # ax.zaxis.set_major_formatter('{x.1e}')
        # Add a color bar which maps values to colors.
        # fig.colorbar(surf, shrink=0.5, aspect=5)
        ax.set(
            xlabel='Fibre cost per km',
            ylabel='Initial fibre cost per km',
            zlabel='z_min (bits/s)',
        )
        ax.view_init(20,30)
        plt.savefig("3d_plot_fibres.png")
        plt.show()

    if data_storage_location_keep_each_loop_detect_init_fibre != None:
        if os.path.isfile(data_storage_location_keep_each_loop_detect_init_fibre + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop_detect_init_fibre + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_detector_cost = last_row_explored["Cost_Detector"].iloc[0]
            current_init_fibre_cost = last_row_explored["C_Init_Fibre_conn_per_km"].iloc[0]
            current_key = last_row_explored["Graph key"].iloc[0]
        else:
            current_key = None
            current_detector_cost = 0.0
            current_init_fibre_cost = 0.0
            dictionary_fieldnames = ["Cost_Detector", "C_Init_Fibre_conn_per_km", "Graph key",
                                     "objective_value"]
            with open(data_storage_location_keep_each_loop_detect_init_fibre + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_key = None
        current_detector_cost = 0.0
        current_init_fibre_cost = 0.0
    # for c_det in np.arange(current_detector_cost, c_detector, 0.1):
    no_solution_list = []
    objective_list = {}

    for c_det in np.arange(current_detector_cost, c_detector, 0.1):
        for c_initial_fibre in np.arange(0.0, c_initial_fibre_conn_per_km, 0.05):
            if np.abs(c_det - current_detector_cost) < 0.0001:
                if c_initial_fibre < current_init_fibre_cost - 0.0001:
                    continue
            for key in g.keys():
                if current_key != None and current_key != key:
                    continue
                elif current_key != None:
                    current_key = None
                    continue
                if key not in no_solution_list:
                    try:
                        prob = cplex.Cplex()
                        optimisation = Trusted_Node_Optimisation_No_Switching(prob, g[key], key_dict[key])
                        sol_dict, prob, time_taken = optimisation.trusted_node_optimisation_program(time_limit=2e2,
                                                                                                    cdevices=c_source + c_det,
                                                                                                    c_dist=0.025,
                                                                                                    c_conn_dist=c_initial_fibre,
                                                                                                    cnetwork=overall_cost,
                                                                                                    Lambda=Lambda,
                                                                                                    Nmax=Lambda,
                                                                                                    Tij=Tij)
                        flow_dict, binary_dict, lambda_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                        objective_value = prob.solution.get_objective_value()
                        ### step 2.
                        # if key == 2:
                        #     plot_graph(g[key], delta_dict=delta_dict, save_extra=f"{c_det}_{c_initial_fibre}")
                        if c_det not in objective_list.keys():
                            objective_list[c_det] = {c_initial_fibre: {key: objective_value}}
                        elif c_initial_fibre not in objective_list[c_det].keys():
                            objective_list[c_det][c_initial_fibre] = {key: objective_value}
                        else:
                            objective_list[c_det][c_initial_fibre][key] = objective_value
                        if data_storage_location_keep_each_loop_detect_init_fibre != None:
                            dictionary = [
                                {"Cost_Detector": c_det, "C_Init_Fibre_conn_per_km": c_initial_fibre,
                                 "Graph key": key, "objective_value": objective_value}]
                            dictionary_fieldnames = ["Cost_Detector", "C_Init_Fibre_conn_per_km", "Graph key",
                                                     "objective_value"]
                            if os.path.isfile(data_storage_location_keep_each_loop_detect_init_fibre + '.csv'):
                                with open(data_storage_location_keep_each_loop_detect_init_fibre + '.csv',
                                          mode='a') as csv_file:
                                    writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                    writer.writerows(dictionary)
                            else:
                                with open(data_storage_location_keep_each_loop_detect_init_fibre + '.csv',
                                          mode='a') as csv_file:
                                    writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                    writer.writeheader()
                                    writer.writerows(dictionary)
                    except:
                        no_solution_list.append(key)
                        continue
    if data_storage_location_keep_each_loop_detect_init_fibre != None:
        plot_information = pd.read_csv(data_storage_location_keep_each_loop_detect_init_fibre + ".csv")
        for index, row in plot_information.iterrows():
            if round(row["Cost_Detector"],1) not in objective_list.keys():
                objective_list[round(row["Cost_Detector"],1)] = {
                    round(row["C_Init_Fibre_conn_per_km"],2): {row["Graph key"]: row["objective_value"]}}
            elif round(row["C_Init_Fibre_conn_per_km"],2)  not in objective_list[round(row["Cost_Detector"],1)].keys():
                objective_list[round(row["Cost_Detector"],1)][round(row["C_Init_Fibre_conn_per_km"],2)] = {
                    row["Graph key"]: row["objective_value"]}
            else:
                objective_list[round(row["Cost_Detector"],1)][round(row["C_Init_Fibre_conn_per_km"],2)][row["Graph key"]] = row[
                    "objective_value"]
        #### figure out what investigation to do here
        #### Firstly- this is relative to the cost of the source... c_{det}, c_{fibre_per_km}, c_{init_fibre}
        #### 1. Look at how each of these terms variation affects the overall graph cost - which one has the biggest effect
        #### 2. Investigate relationships between the cost terms, e.g. if the fibre placement cost scales the optimal
        #### solution will inevitably result in different graphs. How do the structures of these graphs vary with this cost
        #### And same for detector cost
        #### 3. Plot a set of 3-D graphs of max min capacity based on 2 parameters...

        #### Standard price for source: 1.5k, Standard price for detectors 120k, per fibre km: 3.8k, init cost to bury: 32.7k
        #### if we set the price of source to 0.01 --> det price: 0.8, per fibre km: 0.025, init cost to bury: 0.22 (set this to 0.2)

        #### Program for 1.
        # objectives_detector_cost = {}
        # objective_costs_to_normalise  = {}
        # for cost in objective_list.keys():
        #     for key in objective_list[cost][0.025][0.2].keys():
        #         if cost < 0.00001:
        #             objective_costs_to_normalise[key] = objective_list[cost][0.025][0.2][key]
        #         if objective_costs_to_normalise[key] > 0.00001:
        #             if cost not in objectives_detector_cost.keys():
        #                 objectives_detector_cost[cost] = [objective_list[cost][0.025][0.2][key]/objective_costs_to_normalise[key]]
        #             else:
        #                 objectives_detector_cost[cost].extend([objective_list[cost][0.025][0.2][key]/objective_costs_to_normalise[key]])
        # mean_values = []
        # std_values = []
        # x = []
        # for cost in objective_costs_to_normalise.keys():
        #     mean_values.append(np.mean(objectives_detector_cost[cost]))
        #     std_values.append(np.std(objectives_detector_cost[cost]))
        #     x.append(cost)
        #
        # plt.errorbar(x, mean_values, yerr=std_values, color="r", capsize=10)
        # plt.xlabel("Cost of Detectors", fontsize=10)
        # plt.ylabel("z_min (bits/s)", fontsize=10)
        # plt.savefig("detector_cost_investigation.png")
        # plt.show()


        objectives_detector = {}
        objective_costs_to_normalise = {}
        for cost in objective_list.keys():
            for key in objective_list[cost][0.2].keys():
                if cost < 0.00001:
                    objective_costs_to_normalise[key] = objective_list[cost][0.2][key]
                if objective_costs_to_normalise[key] > 0.00001:
                    if cost not in objectives_detector.keys():
                        objectives_detector[cost] = [
                            objective_list[cost][0.2][key] / objective_costs_to_normalise[key]]
                    else:
                        objectives_detector[cost].extend(
                            [objective_list[cost][0.2][key] / objective_costs_to_normalise[key]])
        mean_values = []
        std_values = []
        x = []
        for cost in objectives_detector.keys():
            mean_values.append(np.mean(objectives_detector[cost]))
            std_values.append(np.std(objectives_detector[cost]))
            x.append(cost)

        plt.errorbar(x, mean_values, yerr=std_values, color="r", capsize=10)
        plt.xlabel("Detector Cost", fontsize=10)
        plt.ylabel("z_min  / z_min at zero cost", fontsize=10)
        plt.savefig("detector_cost_investigation.png")
        plt.show()


        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        Y = np.arange(0.0, c_initial_fibre_conn_per_km, 0.05)
        X = np.arange(0.0, c_detector, 0.1)
        X, Y = np.meshgrid(X, Y)
        Z = np.empty((len(X), len(X[0])))
        for i in range(len(X)):
            for j in range(len(X[i])):
                Z[i][j] = np.mean(list(objective_list[round(X[i][j], 3)][round(Y[i][j], 2)].values()))
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        #
        zmin, zmax = Z.min(), Z.max()
        sig_fig_min = int(np.floor(np.log10(abs(zmin))))
        sig_fig_max = int(np.floor(np.log10(abs(zmax))))
        ax.set_zlim(np.floor(zmin / (10 ** sig_fig_min)) * (10 ** sig_fig_min),
                    np.ceil(zmax / (10 ** (sig_fig_max - 1))) * (10 ** (sig_fig_max - 1)))
        ax.zaxis.set_major_locator(LinearLocator(10))
        # A StrMethodFormatter is used automatically
        # ax.zaxis.set_major_formatter('{x:2e}')
        # Add a color bar which maps values to colors.
        # fig.colorbar(surf, shrink=0.5, aspect=5)
        ax.set(
            xlabel='Detector cost',
            ylabel='Initial fibre cost per km',
            zlabel='z_min (bits/s)',
        )
        ax.view_init(20, 30)
        plt.savefig("3d_plot_init_fibres_detectors.png")
        plt.show()

        if data_storage_location_keep_each_loop_detect_fibre != None:
            if os.path.isfile(data_storage_location_keep_each_loop_detect_fibre + '.csv'):
                plot_information = pd.read_csv(data_storage_location_keep_each_loop_detect_fibre + ".csv")
                last_row_explored = plot_information.iloc[[-1]]
                current_detector_cost = last_row_explored["Cost_Detector"].iloc[0]
                current_fibre_cost = last_row_explored["Cost_fibre_per_km"].iloc[0]
                current_key = last_row_explored["Graph key"].iloc[0]
            else:
                current_key = None
                current_detector_cost = 0.0
                current_fibre_cost = 0.0
                dictionary_fieldnames = ["Cost_Detector", "Cost_fibre_per_km", "Graph key",
                                         "objective_value"]
                with open(data_storage_location_keep_each_loop_detect_fibre + '.csv', mode='a') as csv_file:
                    writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                    writer.writeheader()
        else:
            current_key = None
            current_detector_cost = 0.0
            current_fibre_cost = 0.0
        # for c_det in np.arange(current_detector_cost, c_detector, 0.1):
        no_solution_list = []
        objective_list = {}

        for c_det in np.arange(current_detector_cost, c_detector, 0.1):
            for c_fibre in np.arange(0.0, c_fibre_per_km, 0.025):
                if np.abs(c_det - current_detector_cost) < 0.0001:
                    if c_fibre < current_fibre_cost - 0.0001:
                        continue

                for key in g.keys():
                    if current_key != None and current_key != key:
                        continue
                    elif current_key != None:
                        current_key = None
                        continue
                    if key not in no_solution_list:
                        try:
                            prob = cplex.Cplex()
                            optimisation = Trusted_Node_Optimisation_No_Switching(prob, g[key], key_dict[key])
                            sol_dict, prob, time_taken = optimisation.trusted_node_optimisation_program(time_limit=2e2,
                                                                                                        cdevices=c_source + c_det,
                                                                                                        c_dist=c_fibre,
                                                                                                        c_conn_dist=0.2,
                                                                                                        cnetwork=overall_cost,
                                                                                                        Lambda=Lambda,
                                                                                                        Nmax=Lambda,
                                                                                                        Tij=Tij)
                            flow_dict, binary_dict, lambda_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                            objective_value = prob.solution.get_objective_value()
                            ### step 2.
                            # if key == 2:
                            #     plot_graph(g[key], delta_dict=delta_dict, save_extra=f"{c_det}_{c_fibre}")
                            if c_det not in objective_list.keys():
                                objective_list[c_det] = {c_fibre: {key: objective_value}}
                            elif c_fibre not in objective_list[c_det].keys():
                                objective_list[c_det][c_fibre] = {key: objective_value}
                            else:
                                objective_list[c_det][c_fibre][key] = objective_value
                            if data_storage_location_keep_each_loop_detect_fibre != None:
                                dictionary = [
                                    {"Cost_Detector": c_det, "Cost_fibre_per_km": c_fibre,
                                     "Graph key": key, "objective_value": objective_value}]
                                dictionary_fieldnames = ["Cost_Detector", "Cost_fibre_per_km", "Graph key",
                                                         "objective_value"]
                                if os.path.isfile(data_storage_location_keep_each_loop_detect_fibre + '.csv'):
                                    with open(data_storage_location_keep_each_loop_detect_fibre + '.csv',
                                              mode='a') as csv_file:
                                        writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                        writer.writerows(dictionary)
                                else:
                                    with open(data_storage_location_keep_each_loop_detect_fibre + '.csv',
                                              mode='a') as csv_file:
                                        writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                        writer.writeheader()
                                        writer.writerows(dictionary)
                        except:
                            no_solution_list.append(key)
                            continue
        if data_storage_location_keep_each_loop_detect_fibre != None:
            plot_information = pd.read_csv(data_storage_location_keep_each_loop_detect_fibre + ".csv")
            for index, row in plot_information.iterrows():
                if round(row["Cost_Detector"],1) not in objective_list.keys():
                    objective_list[round(row["Cost_Detector"],1)] = {
                        round(row["Cost_fibre_per_km"],3): {row["Graph key"]: row["objective_value"]}}
                elif round(row["Cost_fibre_per_km"],3)  not in objective_list[round(row["Cost_Detector"],1)].keys():
                    objective_list[round(row["Cost_Detector"],1)][round(row["Cost_fibre_per_km"],3)] = {
                        row["Graph key"]: row["objective_value"]}
                else:
                    objective_list[round(row["Cost_Detector"],1)][round(row["Cost_fibre_per_km"],3)][row["Graph key"]] = row[
                        "objective_value"]


        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        Y = np.arange(0.0, 0.2, 0.025)
        X = np.arange(0.0, c_detector, 0.1)
        X, Y = np.meshgrid(X, Y)
        Z = np.empty((len(X), len(X[0])))
        for i in range(len(X)):
            for j in range(len(X[i])):
                Z[i][j] = np.mean(list(objective_list[round(X[i][j], 1)][round(Y[i][j], 3)].values()))
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        #
        zmin, zmax = Z.min(), Z.max()
        sig_fig_min = int(np.floor(np.log10(abs(zmin))))
        sig_fig_max = int(np.floor(np.log10(abs(zmax))))
        ax.set_zlim(np.floor(zmin / (10 ** sig_fig_min)) * (10 ** sig_fig_min),
                    np.ceil(zmax / (10 ** (sig_fig_max - 1))) * (10 ** (sig_fig_max - 1)))
        ax.zaxis.set_major_locator(LinearLocator(10))
        # A StrMethodFormatter is used automatically
        # ax.zaxis.set_major_formatter('{x:.02f}')
        # Add a color bar which maps values to colors.
        # fig.colorbar(surf, shrink=0.5, aspect=5)
        ax.set(
            xlabel='Detector cost',
            ylabel='Fibre cost per km',
            zlabel='z_min (bits/s)',
        )
        ax.view_init(20, 30)
        plt.savefig("3d_plot_fibres_detectors.png")
        plt.show()


def cost_detector_investigation_switching(c_source, c_detector, c_fibre_per_km, c_initial_fibre_conn_per_km, overall_cost,frac_switch, Lambda, Tij, node_file_path, edge_file_path, edge_no_switching_path, key_dict_file_path, data_storage_location_keep_each_loop_no_switching,data_storage_location_keep_each_loop_switching):
    key_dict, g = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path,
                                                                    edge_file_path=edge_file_path,
                                                                    key_dict_file_path=key_dict_file_path)
    if data_storage_location_keep_each_loop_switching != None:
        if os.path.isfile(data_storage_location_keep_each_loop_switching + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop_switching + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_detector = last_row_explored["C_detector"].iloc[0]
            current_key = last_row_explored["Graph key"].iloc[0]
        else:
            current_key = None
            current_detector = 0.0
            dictionary_fieldnames = ["C_detector", "Graph key", "objective_value"]
            with open(data_storage_location_keep_each_loop_switching + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_key = None
        current_detector = 0.0
    # for c_det in np.arange(current_detector_cost, c_detector, 0.1):
    no_solution_list = []
    objective_list = {}

    for c_det in np.arange(current_detector, c_detector, 0.1):
        for key in g.keys():
            if current_key != None and current_key != key:
                continue
            elif current_key != None:
                current_key = None
                continue
            if key not in no_solution_list and key < 10:
                try:
                    key_dict_new = switched_fully_trusted_network.make_key_dict_bidirectional(key_dict[key])
                    prob = cplex.Cplex()
                    optimisation = switched_fully_trusted_network.Switched_Trusted_Node_Optimisation_Constant_Switching_Ratio(
                        prob, g[key], key_dict_new)
                    sol_dict, prob = optimisation.trusted_node_optimisation_program(time_limit=2e2, csource=c_source,
                                                                                    cdet=c_det,
                                                                                    c_dist=c_fibre_per_km,
                                                                                    c_conn_dist=c_initial_fibre_conn_per_km,
                                                                                    Lambda=Lambda, cgraph=overall_cost,
                                                                                    f_switch=frac_switch, Nmax=Lambda,
                                                                                    Tij=Tij, no_channels=1)
                    flow_dict, numbers_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                    objective_value = prob.solution.get_objective_value()
                    ### step 2.
                    # if key == 2:
                    #     plot_graph(g[key], delta_dict=delta_dict, save_extra=f"{c_fibre}_{c_initial_fibre}")
                    if c_det not in objective_list.keys():
                        objective_list[c_det] =  {key: objective_value}
                    else:
                        objective_list[c_det][key] = objective_value
                    if data_storage_location_keep_each_loop_switching != None:
                        dictionary = [
                            {"C_detector": c_det, "Graph key": key, "objective_value": objective_value}]
                        dictionary_fieldnames = ["C_detector", "Graph key", "objective_value"]
                        if os.path.isfile(data_storage_location_keep_each_loop_switching + '.csv'):
                            with open(data_storage_location_keep_each_loop_switching + '.csv', mode='a') as csv_file:
                                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                writer.writerows(dictionary)
                        else:
                            with open(data_storage_location_keep_each_loop_switching + '.csv', mode='a') as csv_file:
                                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                writer.writeheader()
                                writer.writerows(dictionary)
                except:
                    no_solution_list.append(key)
                    continue
    if data_storage_location_keep_each_loop_switching != None:
        plot_information = pd.read_csv(data_storage_location_keep_each_loop_switching + ".csv")
        for index, row in plot_information.iterrows():
            if round(row["C_detector"],1) not in objective_list.keys():
                objective_list[round(row["C_detector"],1)] = {row["Graph key"]: row["objective_value"]}
            else:
                objective_list[round(row["C_detector"],1)][row["Graph key"]] = row["objective_value"]
    objective_no_switching_list = {}

    key_dict_no_switch, g_no_switch = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path,
                                                                    edge_file_path=edge_no_switching_path,
                                                                    key_dict_file_path=key_dict_file_path)

    if data_storage_location_keep_each_loop_no_switching != None:
        if os.path.isfile(data_storage_location_keep_each_loop_no_switching + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop_no_switching + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_detector = last_row_explored["C_detector"].iloc[0]
            current_key = last_row_explored["Graph key"].iloc[0]
        else:
            current_key = None
            current_detector = 0.0
            dictionary_fieldnames = ["C_detector", "Graph key", "objective_value"]
            with open(data_storage_location_keep_each_loop_no_switching + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_key = None
        current_detector = 0.0
    # for c_det in np.arange(current_detector_cost, c_detector, 0.1):
    no_solution_list = []

    for c_det in np.arange(current_detector, c_detector, 0.1):
        for key in g_no_switch.keys():
            if current_key != None and current_key != key:
                continue
            elif current_key != None:
                current_key = None
                continue
            if key not in no_solution_list and key < 10:
                try:
                    prob = cplex.Cplex()
                    optimisation = Trusted_Node_Optimisation_No_Switching(prob, g_no_switch[key], key_dict_no_switch[key])
                    sol_dict, prob, time_taken = optimisation.trusted_node_optimisation_program(time_limit=2e2,
                                                                                                cdevices=c_source + c_det,
                                                                                                c_dist=c_fibre_per_km,
                                                                                                c_conn_dist=c_initial_fibre_conn_per_km,
                                                                                                cnetwork=overall_cost,
                                                                                                Lambda=Lambda,
                                                                                                Nmax=Lambda, Tij=Tij,
                                                                                                no_channels=1)
                    flow_dict, binary_dict, lambda_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                    objective_value = prob.solution.get_objective_value()
                    ### step 2.
                    # if key == 2:
                    #     plot_graph(g[key], delta_dict=delta_dict, save_extra=f"{c_fibre}_{c_initial_fibre}")
                    if c_det not in objective_list.keys():
                        objective_list[c_det] =  {key: objective_value}
                    else:
                        objective_list[c_det][key] = objective_value
                    if data_storage_location_keep_each_loop_no_switching != None:
                        dictionary = [
                            {"C_detector": c_det, "Graph key": key, "objective_value": objective_value}]
                        dictionary_fieldnames = ["C_detector", "Graph key", "objective_value"]
                        if os.path.isfile(data_storage_location_keep_each_loop_no_switching + '.csv'):
                            with open(data_storage_location_keep_each_loop_no_switching + '.csv', mode='a') as csv_file:
                                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                writer.writerows(dictionary)
                        else:
                            with open(data_storage_location_keep_each_loop_no_switching + '.csv', mode='a') as csv_file:
                                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                writer.writeheader()
                                writer.writerows(dictionary)
                except:
                    no_solution_list.append(key)
                    continue

    if data_storage_location_keep_each_loop_no_switching != None:
        plot_information = pd.read_csv(data_storage_location_keep_each_loop_no_switching + ".csv")
        for index, row in plot_information.iterrows():
            if round(row["C_detector"], 1) not in objective_no_switching_list.keys():
                objective_no_switching_list[round(row["C_detector"], 1)] = {row["Graph key"]: row["objective_value"]}
            else:
                objective_no_switching_list[round(row["C_detector"], 1)][row["Graph key"]] = row["objective_value"]



    objectives_fractional_detector_differences = {}
    for det_cost in objective_list.keys():
        if det_cost < 1.55:
            for key in objective_list[det_cost].keys():
                if det_cost not in objectives_fractional_detector_differences.keys():
                    objectives_fractional_detector_differences[det_cost] = [
                        (objective_list[det_cost][key]) / objective_no_switching_list[det_cost][key] - 1]
                else:
                    objectives_fractional_detector_differences[det_cost].extend(
                        [ (objective_list[det_cost][key]) / objective_no_switching_list[det_cost][key] - 1])
    mean_values = []
    std_values = []
    x = []
    for cost in objectives_fractional_detector_differences.keys():
        mean_values.append(np.mean(objectives_fractional_detector_differences[cost]))
        std_values.append(np.std(objectives_fractional_detector_differences[cost]))
        x.append(cost)

    plt.errorbar(x, mean_values, yerr=std_values, color="r", capsize=10)
    plt.xlabel("Cost of Detectors", fontsize=10)
    plt.ylabel("(z_min(switch) - z_min(No Switch)) / z_min(No Switch)", fontsize=10)
    plt.savefig("detector_cost_switching_comaprison_0_cost_init_fibre_20_Cost.png")
    plt.show()


def cost_overall_investigation_switching(c_source, c_detector, c_fibre_per_km, c_initial_fibre_conn_per_km,
                                          overall_cost, frac_switch, Lambda, Tij, node_file_path, edge_file_path,
                                          edge_no_switching_path, key_dict_file_path,
                                          data_storage_location_keep_each_loop_no_switching,
                                          data_storage_location_keep_each_loop_switching):
    key_dict, g = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path,
                                                                    edge_file_path=edge_file_path,
                                                                    key_dict_file_path=key_dict_file_path)
    if data_storage_location_keep_each_loop_switching != None:
        if os.path.isfile(data_storage_location_keep_each_loop_switching + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop_switching + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_total_cost = last_row_explored["C_tot"].iloc[0]
            current_key = last_row_explored["Graph key"].iloc[0]
        else:
            current_key = None
            current_total_cost = 100.0
            dictionary_fieldnames = ["C_tot", "Graph key", "objective_value"]
            with open(data_storage_location_keep_each_loop_switching + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_key = None
        current_total_cost = 100
    # for c_det in np.arange(current_detector_cost, c_detector, 0.1):
    no_solution_list = []
    objective_list = {}
    for c_overall in np.arange(current_total_cost, overall_cost, 10):
        for key in g.keys():
            if current_key != None and current_key != key:
                continue
            elif current_key != None:
                current_key = None
                continue
            try:
                key_dict_new = switched_fully_trusted_network.make_key_dict_bidirectional(key_dict[key])
                prob = cplex.Cplex()
                optimisation = switched_fully_trusted_network.Switched_Trusted_Node_Optimisation_Constant_Switching_Ratio(
                    prob, g[key], key_dict_new)
                sol_dict, prob = optimisation.trusted_node_optimisation_program(time_limit=2e2, csource=c_source,
                                                                                cdet=c_detector,
                                                                                c_dist=c_fibre_per_km,
                                                                                c_conn_dist=c_initial_fibre_conn_per_km,
                                                                                Lambda=Lambda, cgraph=c_overall,
                                                                                f_switch=frac_switch, Nmax=Lambda,
                                                                                Tij=Tij, no_channels=1)
                flow_dict, numbers_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                objective_value = prob.solution.get_objective_value()
                ### step 2.
                # if key == 2:
                #     plot_graph(g[key], delta_dict=delta_dict, save_extra=f"{c_fibre}_{c_initial_fibre}")
                if c_overall not in objective_list.keys():
                    objective_list[c_overall] = {key: objective_value}
                else:
                    objective_list[c_overall][key] = objective_value
                if data_storage_location_keep_each_loop_switching != None:
                    dictionary = [
                        {"C_tot": c_overall, "Graph key": key, "objective_value": objective_value}]
                    dictionary_fieldnames = ["C_tot", "Graph key", "objective_value"]
                    if os.path.isfile(data_storage_location_keep_each_loop_switching + '.csv'):
                        with open(data_storage_location_keep_each_loop_switching + '.csv', mode='a') as csv_file:
                            writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                            writer.writerows(dictionary)
                    else:
                        with open(data_storage_location_keep_each_loop_switching + '.csv', mode='a') as csv_file:
                            writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                            writer.writeheader()
                            writer.writerows(dictionary)
            except:
                continue
    if data_storage_location_keep_each_loop_switching != None:
        plot_information = pd.read_csv(data_storage_location_keep_each_loop_switching + ".csv")
        for index, row in plot_information.iterrows():
            if round(row["C_tot"], 0) not in objective_list.keys():
                objective_list[round(row["C_tot"], 0)] = {row["Graph key"]: row["objective_value"]}
            else:
                objective_list[round(row["C_tot"], 0)][row["Graph key"]] = row["objective_value"]
    objective_no_switching_list = {}

    key_dict_no_switch, g_no_switch = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path,
                                                                                        edge_file_path=edge_no_switching_path,
                                                                                        key_dict_file_path=key_dict_file_path)

    if data_storage_location_keep_each_loop_no_switching != None:
        if os.path.isfile(data_storage_location_keep_each_loop_no_switching + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop_no_switching + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_total_cost = last_row_explored["C_tot"].iloc[0]
            current_key = last_row_explored["Graph key"].iloc[0]
        else:
            current_key = None
            current_total_cost = 100.0
            dictionary_fieldnames = ["C_tot", "Graph key", "objective_value"]
            with open(data_storage_location_keep_each_loop_no_switching + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_key = None
        current_total_cost = 100.0
    # for c_det in np.arange(current_detector_cost, c_detector, 0.1):
    no_solution_list = []

    for c_overall in np.arange(current_total_cost, overall_cost, 10):
        for key in g_no_switch.keys():
            if current_key != None and current_key != key:
                continue
            elif current_key != None:
                current_key = None
                continue
            try:
                prob = cplex.Cplex()
                optimisation = Trusted_Node_Optimisation_No_Switching(prob, g_no_switch[key],
                                                                      key_dict_no_switch[key])
                sol_dict, prob, time_taken = optimisation.trusted_node_optimisation_program(time_limit=2e2,
                                                                                            cdevices=c_source + c_detector,
                                                                                            c_dist=c_fibre_per_km,
                                                                                            c_conn_dist=c_initial_fibre_conn_per_km,
                                                                                            cnetwork=c_overall,
                                                                                            Lambda=Lambda,
                                                                                            Nmax=Lambda, Tij=Tij,
                                                                                            no_channels=1)
                flow_dict, binary_dict, lambda_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                objective_value = prob.solution.get_objective_value()
                ### step 2.
                # if key == 2:
                #     plot_graph(g[key], delta_dict=delta_dict, save_extra=f"{c_fibre}_{c_initial_fibre}")
                if c_overall not in objective_list.keys():
                    objective_list[c_overall] = {key: objective_value}
                else:
                    objective_list[c_overall][key] = objective_value
                if data_storage_location_keep_each_loop_no_switching != None:
                    dictionary = [
                        {"C_tot": c_overall, "Graph key": key, "objective_value": objective_value}]
                    dictionary_fieldnames = ["C_tot", "Graph key", "objective_value"]
                    if os.path.isfile(data_storage_location_keep_each_loop_no_switching + '.csv'):
                        with open(data_storage_location_keep_each_loop_no_switching + '.csv', mode='a') as csv_file:
                            writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                            writer.writerows(dictionary)
                    else:
                        with open(data_storage_location_keep_each_loop_no_switching + '.csv', mode='a') as csv_file:
                            writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                            writer.writeheader()
                            writer.writerows(dictionary)
            except:
                continue

    if data_storage_location_keep_each_loop_no_switching != None:
        plot_information = pd.read_csv(data_storage_location_keep_each_loop_no_switching + ".csv")
        for index, row in plot_information.iterrows():
            if round(row["C_tot"], 0) not in objective_no_switching_list.keys():
                objective_no_switching_list[round(row["C_tot"], 0)] = {row["Graph key"]: row["objective_value"]}
            else:
                objective_no_switching_list[round(row["C_tot"], 0)][row["Graph key"]] = row["objective_value"]

    objectives_fractional_detector_differences = {}
    for overall_cost in objective_list.keys():
        for key in objective_list[overall_cost].keys():
            if overall_cost not in objectives_fractional_detector_differences.keys():
                objectives_fractional_detector_differences[overall_cost] = [
                    (objective_list[overall_cost][key]) / objective_no_switching_list[overall_cost][key]]
            else:
                objectives_fractional_detector_differences[overall_cost].extend(
                    [(objective_list[overall_cost][key]) / objective_no_switching_list[overall_cost][key]])
    mean_values = []
    std_values = []
    x = []
    for cost in objectives_fractional_detector_differences.keys():
        mean_values.append(np.mean(objectives_fractional_detector_differences[cost]))
        std_values.append(np.std(objectives_fractional_detector_differences[cost]))
        x.append(cost)

    plt.errorbar(x, mean_values, yerr=std_values, color="r", capsize=10)
    plt.xlabel("Overall Network Cost", fontsize=10)
    plt.ylabel("z_min with Switching / z_min without switching", fontsize=10)
    plt.savefig("overall_cost_switching_comaprison.png")
    plt.show()
















def variation_of_fibre_cost_effects_with_number_channels(c_source, c_detector, c_fibre_per_km, c_initial_fibre_conn_per_km, overall_cost, Lambda, Tij, node_file_path, edge_file_path, key_dict_file_path, data_storage_location_keep_each_loop_fibre_init_fibre,data_storage_location_keep_each_loop_init_fibre ):
    key_dict, g = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path,
                                                                    edge_file_path=edge_file_path,
                                                                    key_dict_file_path=key_dict_file_path)
    if data_storage_location_keep_each_loop_fibre_init_fibre != None:
        if os.path.isfile(data_storage_location_keep_each_loop_fibre_init_fibre + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop_fibre_init_fibre + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_fibre_per_km = last_row_explored["Cost_fibre_per_km"].iloc[0]
            current_no_connections = last_row_explored["No_connections"].iloc[0]
            current_key = last_row_explored["Graph key"].iloc[0]
        else:
            current_key = None
            current_fibre_per_km = 0.0
            current_no_connections = 1
            dictionary_fieldnames = ["Cost_fibre_per_km", "No_connections", "Graph key", "objective_value"]
            with open(data_storage_location_keep_each_loop_fibre_init_fibre + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_key = None
        current_fibre_per_km = 0.0
        current_no_connections = 1


    # if data_storage_location_keep_each_loop_init_fibre != None:
    #     if os.path.isfile(data_storage_location_keep_each_loop_init_fibre + '.csv'):
    #         plot_information = pd.read_csv(data_storage_location_keep_each_loop_init_fibre + ".csv")
    #         last_row_explored = plot_information.iloc[-1]
    #         current_init_fibre_per_km = last_row_explored["Cost_inital_fibre_per_km"].iloc[0]
    #         current_key_init = last_row_explored["Graph key"].iloc[0]
    #     else:
    #         current_key_init = None
    #         current_init_fibre_per_km = 0.0
    #         dictionary_fieldnames = ["Cost_inital_fibre_per_km", "Graph key", "objective_value"]
    #         with open(data_storage_location_keep_each_loop_init_fibre + '.csv', mode='a') as csv_file:
    #             writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
    #             writer.writeheader()
    # else:
    #     current_key_init = None
    #     current_init_fibre_per_km = 0.0

    # for c_det in np.arange(current_detector_cost, c_detector, 0.1):
    no_solution_list = []
    objective_list = {}
    for conn_no in range(current_no_connections, 10):
        for c_fibre in np.arange(0.0, c_fibre_per_km, 0.025):
            if c_fibre < current_fibre_per_km - 0.0001 and conn_no == current_no_connections:
                continue
            for key in g.keys():
                if current_key != None and current_key != key and conn_no == current_no_connections:
                    continue
                elif current_key != None:
                    current_key = None
                    continue
                if key not in no_solution_list:
                    try:
                        prob = cplex.Cplex()
                        optimisation = Trusted_Node_Optimisation_No_Switching(prob, g[key], key_dict[key])
                        sol_dict, prob, time_taken = optimisation.trusted_node_optimisation_program(time_limit = 2e2, cdevices = c_source + c_detector, c_dist = c_fibre, c_conn_dist = c_initial_fibre_conn_per_km, cnetwork= overall_cost, Lambda = Lambda, Nmax = Lambda, Tij = Tij,no_channels = conn_no)
                        flow_dict, binary_dict, lambda_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                        objective_value = prob.solution.get_objective_value()
                        ### step 2.
                        # if key == 2:
                        #     plot_graph(g[key], delta_dict=delta_dict, save_extra=f"{c_fibre}_{c_initial_fibre}")
                        if conn_no not in objective_list.keys():
                            objective_list[conn_no] = {c_fibre: {key: objective_value}}
                        elif c_fibre not in objective_list[conn_no].keys():
                            objective_list[conn_no][c_fibre] = {key: objective_value}
                        else:
                            objective_list[conn_no][c_fibre][key] = objective_value
                        if data_storage_location_keep_each_loop_fibre_init_fibre != None:
                            dictionary = [
                                {"Cost_fibre_per_km": c_fibre, "No_connections": conn_no, "Graph key": key, "objective_value": objective_value}]
                            dictionary_fieldnames = ["Cost_fibre_per_km", "No_connections", "Graph key", "objective_value"]
                            if os.path.isfile(data_storage_location_keep_each_loop_fibre_init_fibre + '.csv'):
                                with open(data_storage_location_keep_each_loop_fibre_init_fibre + '.csv', mode='a') as csv_file:
                                    writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                    writer.writerows(dictionary)
                            else:
                                with open(data_storage_location_keep_each_loop_fibre_init_fibre + '.csv', mode='a') as csv_file:
                                    writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                    writer.writeheader()
                                    writer.writerows(dictionary)
                    except:
                        no_solution_list.append(key)
                        continue
        # if conn_no == 1:


    if data_storage_location_keep_each_loop_fibre_init_fibre != None:
        plot_information = pd.read_csv(data_storage_location_keep_each_loop_fibre_init_fibre + ".csv")
        for index, row in plot_information.iterrows():
            if row["No_connections"] not in objective_list.keys():
                objective_list[row["No_connections"]] = {round(row["Cost_fibre_per_km"],3): {row["Graph key"]: row["objective_value"]}}
            elif round(row["Cost_fibre_per_km"],2) not in objective_list[row["No_connections"]].keys():
                objective_list[row["No_connections"]][round(row["Cost_fibre_per_km"],3)] = {row["Graph key"]: row["objective_value"]}
            else:
                objective_list[row["No_connections"]][round(row["Cost_fibre_per_km"],3)][row["Graph key"]] = row["objective_value"]
        #### figure out what investigation to do here
        #### Firstly- this is relative to the cost of the source... c_{det}, c_{fibre_per_km}, c_{init_fibre}
        #### 1. Look at how each of these terms variation affects the overall graph cost - which one has the biggest effect
        #### 2. Investigate relationships between the cost terms, e.g. if the fibre placement cost scales the optimal
        #### solution will inevitably result in different graphs. How do the structures of these graphs vary with this cost
        #### And same for detector cost
        #### 3. Plot a set of 3-D graphs of max min capacity based on 2 parameters...

        #### Standard price for source: 1.5k, Standard price for detectors 120k, per fibre km: 3.8k, init cost to bury: 32.7k
        #### if we set the price of source to 0.01 --> det price: 0.8, per fibre km: 0.025, init cost to bury: 0.22 (set this to 0.2)

        #### Program for 1.
        # objectives_detector_cost = {}
        # objective_costs_to_normalise  = {}
        # for cost in objective_list.keys():
        #     for key in objective_list[cost][0.025][0.2].keys():
        #         if cost < 0.00001:
        #             objective_costs_to_normalise[key] = objective_list[cost][0.025][0.2][key]
        #         if objective_costs_to_normalise[key] > 0.00001:
        #             if cost not in objectives_detector_cost.keys():
        #                 objectives_detector_cost[cost] = [objective_list[cost][0.025][0.2][key]/objective_costs_to_normalise[key]]
        #             else:
        #                 objectives_detector_cost[cost].extend([objective_list[cost][0.025][0.2][key]/objective_costs_to_normalise[key]])
        # mean_values = []
        # std_values = []
        # x = []
        # for cost in objective_costs_to_normalise.keys():
        #     mean_values.append(np.mean(objectives_detector_cost[cost]))
        #     std_values.append(np.std(objectives_detector_cost[cost]))
        #     x.append(cost)
        #
        # plt.errorbar(x, mean_values, yerr=std_values, color="r", capsize=10)
        # plt.xlabel("Cost of Detectors", fontsize=10)
        # plt.ylabel("z_min (bits/s)", fontsize=10)
        # plt.savefig("detector_cost_investigation.png")
        # plt.show()
        for no_connections in objective_list.keys():
            objectives_fibre_cost = {}
            objective_costs_to_normalise = {}
            for cost in objective_list[no_connections].keys():
                for key in objective_list[no_connections][cost].keys():
                    if cost < 0.00001:
                        objective_costs_to_normalise[key] = objective_list[no_connections][cost][key]
                    if objective_costs_to_normalise[key] > 0.00001:
                        if cost not in objectives_fibre_cost.keys():
                            objectives_fibre_cost[cost] = [
                                objective_list[no_connections][cost][key] / objective_costs_to_normalise[key]]
                        else:
                            objectives_fibre_cost[cost].extend(
                                [objective_list[no_connections][cost][key] / objective_costs_to_normalise[key]])
            mean_values = []
            std_values = []
            x = []
            for cost in objectives_fibre_cost.keys():
                mean_values.append(np.mean(objectives_fibre_cost[cost]))
                std_values.append(np.std(objectives_fibre_cost[cost]))
                x.append(cost)

            plt.errorbar(x, mean_values, yerr=std_values, capsize=10, label = f"No. of fibre channels used {no_connections}")
        # ax.set(xlabel = "Cost of Fibre per km", ylabel= "z_min / z_min at zero cost")
        plt.xlabel("Cost of Fibre per km", fontsize=10)
        plt.ylabel("z_min / z_min at zero cost", fontsize=10)
        plt.legend(loc = "best")
        plt.savefig("fibre_cost_investigation.png")
        plt.show()



def f_switch_parameter_sweep(node_file_path, edge_file_path, key_dict_file_path,  node_file_path_switching, edge_file_path_switching, key_dict_file_path_switching,overall_cost,
                         c_source, Tij=None, c_detector=1,c_fibre_per_km=0.1, c_initial_fibre_conn_per_km=0.01, Lambda=100, data_storage_location_keep_each_loop = None, data_storage_location_keep_each_loop_no_switch = None):
    key_dict, g = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path_switching,
                                                                    edge_file_path=edge_file_path_switching,
                                                                    key_dict_file_path=key_dict_file_path_switching)
    key_dict_no_switching, g_no_switching = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path,
                                                                    edge_file_path=edge_file_path,
                                                                    key_dict_file_path=key_dict_file_path)

    # graphs = import_graph_structure(node_information=graph_node_data_file, edge_information=graph_edge_data_file)
    if data_storage_location_keep_each_loop != None:
        if os.path.isfile(data_storage_location_keep_each_loop + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            last_ratio_done = last_row_explored["f_switch"]
            dataframe_of_fswitch_done = plot_information[plot_information["f_switch"] == last_ratio_done.iloc[0]]
            current_key = last_row_explored["Graph key"].iloc[0]
            fswitch_current = last_ratio_done.iloc[0]
        else:
            fswitch_current = None
            current_key = None
            dictionary_fieldnames = ["f_switch", "Graph key", "objective_value"]
            with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        fswitch_current = None
        current_key = None

    no_soln_set = []
    objective_value_at_frac_switch_01 = {}
    objective_values_frac_switch= {}
    for frac_switch in np.arange(start=0.0, stop=0.92, step=0.02):
        if fswitch_current != None:
            if frac_switch != fswitch_current:
                continue
            else:
                fswitch_current = None
        for key in g.keys():
            if current_key != None:
                if current_key == key:
                    current_key = None
                    continue
                else:
                    continue
            if key not in no_soln_set:
                try:
                    key_dict_new = switched_fully_trusted_network.make_key_dict_bidirectional(key_dict[key])
                    prob = cplex.Cplex()
                    optimisation = switched_fully_trusted_network.Switched_Trusted_Node_Optimisation_Constant_Switching_Ratio(
                        prob, g[key], key_dict_new)
                    sol_dict, prob = optimisation.trusted_node_optimisation_program(time_limit=2e2, csource=c_source,
                                                                                    cdet=c_detector,
                                                                                    c_dist=c_fibre_per_km,
                                                                                    c_conn_dist=c_initial_fibre_conn_per_km,
                                                                                    Lambda=Lambda, cgraph=overall_cost,
                                                                                    f_switch=frac_switch, Nmax=Lambda,
                                                                                    Tij=Tij)
                    flow_dict, numbers_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                    objective_value = prob.solution.get_objective_value()
                    if abs(frac_switch - 0.1) < 0.0001:
                        objective_value_at_frac_switch_01[key] = objective_value
                    if frac_switch not in objective_values_frac_switch.keys():
                        objective_values_frac_switch[frac_switch] = [(objective_value, key)]
                    else:
                        objective_values_frac_switch[frac_switch].append((objective_value, key))
                    if data_storage_location_keep_each_loop != None:
                        dictionary = [
                            {"f_switch": frac_switch, "Graph key": key, "objective_value": objective_value}]
                        dictionary_fieldnames = ["f_switch", "Graph key", "objective_value"]
                        if os.path.isfile(data_storage_location_keep_each_loop + '.csv'):
                            with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                writer.writerows(dictionary)
                        else:
                            with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                writer.writeheader()
                                writer.writerows(dictionary)
                except:
                    no_soln_set.append(key)
                    continue

    if data_storage_location_keep_each_loop != None:
        plot_information = pd.read_csv(data_storage_location_keep_each_loop + ".csv")
        for index, row in plot_information.iterrows():
            if float(row["f_switch"]) not in objective_values_frac_switch.keys():
                objective_values_frac_switch[float(row["f_switch"])] = {row["Graph key"]: row["objective_value"]}
            else:
                objective_values_frac_switch[float(row["f_switch"])][row["Graph key"]] = row["objective_value"]
        objective_value_at_frac_switch_01 = {}
        for key in objective_values_frac_switch[0.1]:
            objective_value_at_frac_switch_01[key] = objective_values_frac_switch[0.1][key]

    no_soln_set = []
    if data_storage_location_keep_each_loop_no_switch != None:
        if os.path.isfile(data_storage_location_keep_each_loop_no_switch + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop_no_switch + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_key = last_row_explored["Graph key"].iloc[0]
        else:
            current_key = None
            dictionary_fieldnames = ["Graph key", "objective_value"]
            with open(data_storage_location_keep_each_loop_no_switch + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_key = None

    for key in g_no_switching.keys():
        if current_key != None:
            if current_key == key:
                current_key = None
                continue
            else:
                continue
        if key not in no_soln_set:
            try:
                prob = cplex.Cplex()
                optimisation = Trusted_Node_Optimisation_No_Switching(prob, g[key], key_dict[key])
                sol_dict, prob, time_taken = optimisation.trusted_node_optimisation_program(time_limit=2e2,
                                                                                            cdevices=c_source + c_detector,
                                                                                            c_dist=c_fibre_per_km,
                                                                                            c_conn_dist=c_initial_fibre_conn_per_km,
                                                                                            cnetwork=overall_cost,
                                                                                            Lambda=Lambda,
                                                                                            Nmax=Lambda, Tij=Tij)
                flow_dict, binary_dict, lambda_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                objective_value = prob.solution.get_objective_value()
                if data_storage_location_keep_each_loop_no_switch != None:
                    dictionary = [
                        {"Graph key": key, "objective_value": objective_value}]
                    dictionary_fieldnames = ["Graph key", "objective_value"]
                    if os.path.isfile(data_storage_location_keep_each_loop_no_switch + '.csv'):
                        with open(data_storage_location_keep_each_loop_no_switch + '.csv', mode='a') as csv_file:
                            writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                            writer.writerows(dictionary)
                    else:
                        with open(data_storage_location_keep_each_loop_no_switch + '.csv', mode='a') as csv_file:
                            writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                            writer.writeheader()
                            writer.writerows(dictionary)
            except:
                no_soln_set.append(key)
                continue


    if data_storage_location_keep_each_loop_no_switch != None:
        objective_values_no_switch = {}
        plot_information_no_switching = pd.read_csv(data_storage_location_keep_each_loop_no_switch + ".csv")
        for index, row in plot_information_no_switching.iterrows():
            if row["Graph key"] not in objective_values_no_switch.keys():
                objective_values_no_switch[row["Graph key"]] = [row["objective_value"]]
            else:
                objective_values_no_switch[row["Graph key"]].append(row["objective_value"])
        objective_values = {}
        for frac_switch in objective_values_frac_switch.keys():
            for key in objective_values_frac_switch[frac_switch].keys():
                if frac_switch not in objective_values.keys():
                    objective_values[frac_switch] = [objective_values_frac_switch[frac_switch][key] / objective_values_no_switch[key]]
                else:
                    objective_values[frac_switch].append(objective_values_frac_switch[frac_switch][key] / objective_values_no_switch[key])
        mean_objectives = {}
        std_objectives = {}
        for key in objective_values.keys():
            mean_objectives[key] = np.mean([objective_values[key]])
            std_objectives[key] = np.std([objective_values[key]])
        mean_differences = []
        std_differences = []
        # topologies
        x = []
        for key in mean_objectives.keys():
            mean_differences.append(mean_objectives[key])
            std_differences.append(std_objectives[key])
            x.append(key)
        plt.errorbar(x, mean_differences, yerr=std_differences, color="r", capsize= 0, label = "Max z_min of Network")
        plt.axhline(y=1, color='b', linestyle='-', label = "Max z_min of Network without Switching")
        plt.legend()
        plt.xlabel("Fraction of time calibrating", fontsize=10)
        plt.ylabel("Max z_min with Switching/max z_min without Switching", fontsize=10)
        # plt.legend(loc='upper right', fontsize='medium')
        plt.savefig("l_switch_mesh_topology.png")
        plt.show()
    else:
        objective_values = {}
        for frac_switch in objective_values_frac_switch.keys():
            for key in objective_values_frac_switch[frac_switch]:
                if frac_switch not in objective_values.keys():
                    objective_values[frac_switch] = [objective_values_frac_switch[frac_switch][key] / objective_value_at_frac_switch_01[key]]
                else:
                    objective_values[frac_switch].append(objective_values_frac_switch[frac_switch][key] / objective_value_at_frac_switch_01[key])
        mean_objectives = {}
        std_objectives = {}
        for key in objective_values.keys():
            mean_objectives[key] = np.mean(objective_values[key])
            std_objectives[key] = np.std(objective_values[key])
        mean_differences = []
        std_differences = []
        # topologies
        x = []
        for key in mean_objectives.keys():
            mean_differences.append(mean_objectives[key])
            std_differences.append(std_objectives[key])
            x.append(key)
        plt.errorbar(x, mean_differences, yerr=std_differences, color="r", capsize= 0)
        plt.xlabel("Fraction of time calibrating", fontsize=10)
        plt.ylabel("Cost of Network/Cost of Network at 0.1 Fraction of Time Calibrating" , fontsize=10)
        # plt.legend(loc='upper right', fontsize='medium')
        plt.savefig("frac_switch_mesh_topology_01_comp.png")
        plt.show()



def l_switch_parameter_sweep(node_file_path, edge_file_path, key_dict_file_path,  node_file_path_switching, edge_file_path_switching, key_dict_file_path_switching,overall_cost,
                         c_source, f_switch, Tij=None, c_detector=1,c_fibre_per_km=0.1, c_initial_fibre_conn_per_km=0.01, Lambda=100, data_storage_location_keep_each_loop = None, data_storage_location_keep_each_loop_no_switch = None):

    key_dict_no_switching, g_no_switching = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path,
                                                                    edge_file_path=edge_file_path,
                                                                    key_dict_file_path=key_dict_file_path)

    # graphs = import_graph_structure(node_information=graph_node_data_file, edge_information=graph_edge_data_file)
    if data_storage_location_keep_each_loop != None:
        if os.path.isfile(data_storage_location_keep_each_loop + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            last_ratio_done = last_row_explored["L_switch"]
            dataframe_of_fswitch_done = plot_information[plot_information["L_switch"] == last_ratio_done.iloc[0]]
            current_key = last_row_explored["Graph key"].iloc[0]
            l_current = last_ratio_done.iloc[0]
        else:
            l_current = None
            current_key = None
            dictionary_fieldnames = ["L_switch", "Graph key", "objective_value"]
            with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        l_current = None
        current_key = None

    no_soln_set = []
    objective_value_at_frac_switch_01 = {}
    objective_values_l_switch= {}
    for l_switch in np.arange(start=0.25, stop=3, step=0.25):
        if l_current != None:
            if l_switch != l_current:
                continue
            else:
                l_current = None
        key_dict, g = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path_switching,
                                                                        edge_file_path=edge_file_path_switching + str(int(l_switch*100)) + ".csv",
                                                                        key_dict_file_path=key_dict_file_path_switching)

        for key in g.keys():
            if current_key != None:
                if current_key == key:
                    current_key = None
                    continue
                else:
                    continue
            if key not in no_soln_set:
                try:
                    key_dict_new = switched_fully_trusted_network.make_key_dict_bidirectional(key_dict[key])
                    prob = cplex.Cplex()
                    optimisation = switched_fully_trusted_network.Switched_Trusted_Node_Optimisation_Constant_Switching_Ratio(
                        prob, g[key], key_dict_new)
                    sol_dict, prob = optimisation.trusted_node_optimisation_program(time_limit=2e2, csource=c_source,
                                                                                    cdet=c_detector,
                                                                                    c_dist=c_fibre_per_km,
                                                                                    c_conn_dist=c_initial_fibre_conn_per_km,
                                                                                    Lambda=Lambda, cgraph=overall_cost,
                                                                                    f_switch=f_switch, Nmax=Lambda,
                                                                                    Tij=Tij)
                    flow_dict, numbers_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                    objective_value = prob.solution.get_objective_value()
                    if l_switch not in objective_values_l_switch.keys():
                        objective_values_l_switch[l_switch] ={objective_value: key}
                    else:
                        objective_values_l_switch[l_switch][objective_value] =  key
                    if data_storage_location_keep_each_loop != None:
                        dictionary = [
                            {"L_switch": l_switch, "Graph key": key, "objective_value": objective_value}]
                        dictionary_fieldnames = ["L_switch", "Graph key", "objective_value"]
                        if os.path.isfile(data_storage_location_keep_each_loop + '.csv'):
                            with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                writer.writerows(dictionary)
                        else:
                            with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                writer.writeheader()
                                writer.writerows(dictionary)
                except:
                    no_soln_set.append(key)
                    continue

    if data_storage_location_keep_each_loop != None:
        plot_information = pd.read_csv(data_storage_location_keep_each_loop + ".csv")
        for index, row in plot_information.iterrows():
            if float(row["L_switch"]) not in objective_values_l_switch.keys():
                objective_values_l_switch[float(row["L_switch"])] = {row["Graph key"]: row["objective_value"]}
            else:
                objective_values_l_switch[float(row["L_switch"])][row["Graph key"]] = row["objective_value"]
        objective_value_at_frac_switch_01 = {}

    no_soln_set = []
    if data_storage_location_keep_each_loop_no_switch != None:
        if os.path.isfile(data_storage_location_keep_each_loop_no_switch + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop_no_switch + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_key = last_row_explored["Graph key"].iloc[0]
        else:
            current_key = None
            dictionary_fieldnames = ["Graph key", "objective_value"]
            with open(data_storage_location_keep_each_loop_no_switch + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_key = None

    for key in g_no_switching.keys():
        if current_key != None:
            if current_key == key:
                current_key = None
                continue
            else:
                continue
        if key not in no_soln_set:
            try:
                prob = cplex.Cplex()
                optimisation = Trusted_Node_Optimisation_No_Switching(prob, g[key], key_dict[key])
                sol_dict, prob, time_taken = optimisation.trusted_node_optimisation_program(time_limit=2e2,
                                                                                            cdevices=c_source + c_detector,
                                                                                            c_dist=c_fibre_per_km,
                                                                                            c_conn_dist=c_initial_fibre_conn_per_km,
                                                                                            cnetwork=overall_cost,
                                                                                            Lambda=Lambda,
                                                                                            Nmax=Lambda, Tij=Tij)
                flow_dict, binary_dict, lambda_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                objective_value = prob.solution.get_objective_value()
                if data_storage_location_keep_each_loop_no_switch != None:
                    dictionary = [
                        {"Graph key": key, "objective_value": objective_value}]
                    dictionary_fieldnames = ["Graph key", "objective_value"]
                    if os.path.isfile(data_storage_location_keep_each_loop_no_switch + '.csv'):
                        with open(data_storage_location_keep_each_loop_no_switch + '.csv', mode='a') as csv_file:
                            writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                            writer.writerows(dictionary)
                    else:
                        with open(data_storage_location_keep_each_loop_no_switch + '.csv', mode='a') as csv_file:
                            writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                            writer.writeheader()
                            writer.writerows(dictionary)
            except:
                no_soln_set.append(key)
                continue


    if data_storage_location_keep_each_loop_no_switch != None:
        objective_values_no_switch = {}
        plot_information_no_switching = pd.read_csv(data_storage_location_keep_each_loop_no_switch + ".csv")
        for index, row in plot_information_no_switching.iterrows():
            if row["Graph key"] not in objective_values_no_switch.keys():
                objective_values_no_switch[row["Graph key"]] = [row["objective_value"]]
            else:
                objective_values_no_switch[row["Graph key"]].append(row["objective_value"])
        objective_values = {}
        for l_switch in objective_values_l_switch.keys():
            for key in objective_values_l_switch[l_switch].keys():
                if l_switch not in objective_values.keys():
                    objective_values[l_switch] = [objective_values_l_switch[l_switch][key] / objective_values_no_switch[key]]
                else:
                    objective_values[l_switch].append(objective_values_l_switch[l_switch][key] / objective_values_no_switch[key])
        mean_objectives = {}
        std_objectives = {}
        for key in objective_values.keys():
            mean_objectives[key] = np.mean([objective_values[key]])
            std_objectives[key] = np.std([objective_values[key]])
        mean_differences = []
        std_differences = []
        # topologies
        x = []
        for key in mean_objectives.keys():
            mean_differences.append(mean_objectives[key])
            std_differences.append(std_objectives[key])
            x.append(key)
        plt.errorbar(x, mean_differences, yerr=std_differences, color="r", capsize= 0, label = "Max z_min of Network")
        plt.axhline(y=1, color='b', linestyle='-', label = "Max z_min of Network without Switching")
        plt.legend()
        plt.xlabel("L_switch in dB", fontsize=10)
        plt.ylabel("Max z_min with Switching/max z_min without Switching", fontsize=10)
        # plt.legend(loc='upper right', fontsize='medium')
        plt.savefig("l_switch_mesh_topology.png")
        plt.show()
    else:
        objective_values = {}
        for l_switch in objective_values_l_switch.keys():
            for key in objective_values_l_switch[l_switch]:
                if l_switch not in objective_values.keys():
                    objective_values[l_switch] = [objective_values_l_switch[l_switch][key] / objective_value_at_frac_switch_01[key]]
                else:
                    objective_values[l_switch].append(objective_values_l_switch[l_switch][key] / objective_value_at_frac_switch_01[key])
        mean_objectives = {}
        std_objectives = {}
        for key in objective_values.keys():
            mean_objectives[key] = np.mean(objective_values[key])
            std_objectives[key] = np.std(objective_values[key])
        mean_differences = []
        std_differences = []
        # topologies
        x = []
        for key in mean_objectives.keys():
            mean_differences.append(mean_objectives[key])
            std_differences.append(std_objectives[key])
            x.append(key)
        plt.errorbar(x, mean_differences, yerr=std_differences, color="r", capsize= 0)
        plt.xlabel("L_switch in dB", fontsize=10)
        plt.ylabel("Cost of Network/Cost of Network at 0.1 Fraction of Time Calibrating" , fontsize=10)
        # plt.legend(loc='upper right', fontsize='medium')
        plt.savefig("loss_switch_mesh_topology_01_comp.png")
        plt.show()




def comparison_switch_no_switch_different_no_channels(c_source, c_detector, c_fibre_per_km, c_initial_fibre_conn_per_km, overall_cost_min, overall_cost_max,fswitch, Lambda, Tij, node_file_path, edge_file_path, key_dict_file_path,node_file_path_switching, edge_file_path_switching, key_dict_file_path_switching, data_storage_location_keep_each_loop = None, data_storage_location_keep_each_loop_switching = None):
    key_dict, g = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path,
                                                                    edge_file_path=edge_file_path,
                                                                    key_dict_file_path=key_dict_file_path)
    if data_storage_location_keep_each_loop != None:
        if os.path.isfile(data_storage_location_keep_each_loop + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_key = last_row_explored["Graph key"].iloc[0]
            current_channel = last_row_explored["channel_no"].iloc[0]
            current_cost = last_row_explored["overall_cost"].iloc[0]
        else:
            current_key = None
            current_channel = 1
            current_cost = None
            dictionary_fieldnames = ["Graph key", "channel_no","overall_cost", "objective_value"]
            with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_key = None
        current_channel = 1
        current_cost = None

    objective_values = {}
    no_solution_list = []
    for channel_number in np.arange(current_channel,6,1):
        for overall_cost in np.arange(overall_cost_min, overall_cost_max, 10):
            if current_cost != None and current_cost != overall_cost:
                continue
            elif current_cost != None:
                current_cost = None
            for key in g.keys():
                if current_key != None and key != current_key:
                    continue
                elif current_key != None:
                    current_key = None
                    continue
                if key not in no_solution_list:
                    try:
                        prob = cplex.Cplex()
                        optimisation = Trusted_Node_Optimisation_No_Switching(prob, g[key], key_dict[key])
                        sol_dict, prob, time_taken = optimisation.trusted_node_optimisation_program(time_limit=2e2,
                                                                                        cdevices=c_source + c_detector,
                                                                                        c_dist=c_fibre_per_km,
                                                                                        c_conn_dist=c_initial_fibre_conn_per_km,
                                                                                        cnetwork=overall_cost,
                                                                                        Lambda=Lambda, Nmax = Lambda, Tij=Tij, no_channels = channel_number)
                        flow_dict, binary_dict, lambda_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                        objective_value = prob.solution.get_objective_value()
                        # if key == 2:
                        #     plot_graph(g[key], delta_dict=delta_dict, save_extra=f"comparison_no_switching_{key}_cost_{overall_cost}")
                        if channel_number not in objective_values.keys():
                            objective_values[channel_number]= {overall_cost: {key:objective_value}}
                        elif overall_cost not in objective_values[channel_number].keys():
                            objective_values[channel_number][overall_cost]= {key:objective_value}
                        else:
                            objective_values[channel_number][overall_cost][key] = objective_value

                        if data_storage_location_keep_each_loop != None:
                            dictionary = [
                                {"Graph key": key, "channel_no" : channel_number, "overall_cost": overall_cost,
                                 "objective_value": objective_value}]
                            dictionary_fieldnames = ["Graph key","channel_no","overall_cost", "objective_value"]
                            if os.path.isfile(data_storage_location_keep_each_loop + '.csv'):
                                with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                                    writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                    writer.writerows(dictionary)
                            else:
                                with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                                    writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                    writer.writeheader()
                                    writer.writerows(dictionary)
                    except:
                        continue
    key_dict_switching, g_switching = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path_switching,
                                                                    edge_file_path=edge_file_path_switching,
                                                                    key_dict_file_path=key_dict_file_path_switching)
    if data_storage_location_keep_each_loop_switching != None:
        if os.path.isfile(data_storage_location_keep_each_loop_switching + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop_switching + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_key = last_row_explored["Graph key"].iloc[0]
            current_channel = last_row_explored["channel_no"].iloc[0]
            current_cost = last_row_explored["overall_cost"].iloc[0]
        else:
            current_key = None
            current_channel = 1
            current_cost = None
            dictionary_fieldnames = ["Graph key","channel_no","overall_cost", "objective_value"]
            with open(data_storage_location_keep_each_loop_switching + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_key = None
        current_channel = 1
        current_cost = None

    objective_values_switching = {}
    no_solution_list_switching = []
    for channel_number in np.arange(current_channel, 6, 1):
        for overall_cost in np.arange(overall_cost_min, overall_cost_max, 10):
            if current_cost != None and current_cost != overall_cost:
                continue
            elif current_cost != None:
                current_cost = None
            for key in g_switching.keys():
                if current_key != None and key != current_key:
                    continue
                elif current_key != None:
                    current_key = None
                    continue
                if key not in no_solution_list_switching:
                    try:
                        key_dict_new = switched_fully_trusted_network.make_key_dict_bidirectional(key_dict[key])
                        prob = cplex.Cplex()
                        optimisation = switched_fully_trusted_network.Switched_Trusted_Node_Optimisation_Constant_Switching_Ratio(prob, g[key], key_dict_new)
                        sol_dict, prob = optimisation.trusted_node_optimisation_program(time_limit=2e2, csource=c_source, cdet=c_detector,
                                                                        c_dist=c_fibre_per_km, c_conn_dist=c_initial_fibre_conn_per_km, Lambda=Lambda, cgraph=overall_cost,
                                                                        f_switch=fswitch, Nmax = Lambda, Tij = Tij, no_channels = channel_number)
                        flow_dict, numbers_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                        objective_value = prob.solution.get_objective_value()
                        # if key == 2:
                        #     plot_graph(g[key], delta_dict=delta_dict, save_extra=f"comparison_switching_{key}_cost_{overall_cost}")
                        if channel_number not in objective_values_switching.keys():
                            objective_values_switching[channel_number]= {overall_cost: {key:objective_value}}
                        elif overall_cost not in objective_values_switching[channel_number].keys():
                            objective_values_switching[channel_number][overall_cost]= {key:objective_value}
                        else:
                            objective_values_switching[channel_number][overall_cost][key] = objective_value

                        if data_storage_location_keep_each_loop_switching != None:
                            dictionary = [
                                {"Graph key": key, "channel_no": channel_number, "overall_cost": overall_cost,
                                 "objective_value": objective_value}]
                            dictionary_fieldnames = ["Graph key","channel_no", "overall_cost", "objective_value"]
                            if os.path.isfile(data_storage_location_keep_each_loop_switching + '.csv'):
                                with open(data_storage_location_keep_each_loop_switching + '.csv', mode='a') as csv_file:
                                    writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                    writer.writerows(dictionary)
                            else:
                                with open(data_storage_location_keep_each_loop_switching + '.csv', mode='a') as csv_file:
                                    writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                    writer.writeheader()
                                    writer.writerows(dictionary)
                    except:
                        continue
    if data_storage_location_keep_each_loop != None:
        plot_information = pd.read_csv(data_storage_location_keep_each_loop + ".csv")
        for index, row in plot_information.iterrows():
            if row["channel_no"] not in objective_values.keys():
                objective_values[row["channel_no"]] = {row["overall_cost"] : {row["Graph key"]: row["objective_value"]}}
            if row["overall_cost"] not in objective_values[row["channel_no"]].keys():
                objective_values[row["channel_no"]][row["overall_cost"]] = {row["Graph key"]: row["objective_value"]}
            else:
                objective_values[row["channel_no"]][row["overall_cost"]][row["Graph key"]] = row["objective_value"]
    if data_storage_location_keep_each_loop_switching != None:
        plot_information = pd.read_csv(data_storage_location_keep_each_loop_switching + ".csv")
        for index, row in plot_information.iterrows():
            if row["channel_no"] not in objective_values_switching.keys():
                objective_values_switching[row["channel_no"]] = {row["overall_cost"]: {row["Graph key"]: row["objective_value"]}}
            if row["overall_cost"] not in objective_values_switching[row["channel_no"]].keys():
                objective_values_switching[row["channel_no"]][row["overall_cost"]] = {row["Graph key"]: row["objective_value"]}
            else:
                objective_values_switching[row["channel_no"]][row["overall_cost"]][row["Graph key"]] = row["objective_value"]

    objective_differences = {}
    for channel_no in objective_values.keys():
        for cost in objective_values[channel_no].keys():
            for key in objective_values[channel_no][cost].keys():
                if objective_values[channel_no][cost][key] > 0.0001:
                    if channel_no not in objective_differences.keys():
                        objective_differences[channel_no] ={cost:  [(objective_values_switching[channel_no][cost][key] - objective_values[channel_no][cost][key]) / objective_values[channel_no][cost][key]]}
                    elif cost not in objective_differences[channel_no].keys():
                        objective_differences[channel_no][cost] = [(objective_values_switching[channel_no][cost][key] - objective_values[channel_no][cost][key]) / objective_values[channel_no][cost][key]]
                    else:
                        objective_differences[channel_no][cost].append((objective_values_switching[channel_no][cost][key] - objective_values[channel_no][cost][key])/objective_values[channel_no][cost][key])
    mean_objectives = {}
    std_objectives = {}
    for key in objective_differences.keys():
        for cost in objective_differences[key].keys():
            if key not in mean_objectives.keys():
                mean_objectives[key] = {cost: np.mean(objective_differences[key][cost])}
                std_objectives[key] = {cost: np.std(objective_differences[key][cost])}
            else:
                mean_objectives[key][cost] = np.mean(objective_differences[key][cost])
                std_objectives[key][cost] = np.std(objective_differences[key][cost])
    for no_channels in mean_objectives.keys():
        mean_differences = []
        std_differences = []
        # topologies
        x = []
        for key in mean_objectives[no_channels].keys():
            mean_differences.append(mean_objectives[no_channels][key])
            std_differences.append(std_objectives[no_channels][key])
            x.append(key)
        plt.errorbar(x, mean_differences, yerr=std_differences, label = str(no_channels) + " channels per fibre used")
        plt.xlabel("Overall Cost of the Graph", fontsize=10)
        plt.ylabel("(z_min(switch)-z_min(no switch))/z_min(no switch)", fontsize=10)
    plt.legend(loc='upper right', fontsize='medium')
    plt.savefig("switching_no_switching_comparison.png")
    plt.show()








def monte_carlo_T_i_j(c_source, c_detector, c_fibre_per_km, c_initial_fibre_conn_per_km, overall_cost, Lambda, Tij_var, node_file_path, edge_file_path, key_dict_file_path, num_iterations, data_storage_location_keep_each_loop):
    ### Monte carlo simulation of looking at various different Tij values and determining the optimal throughput for a given graph
    ### wish to look at how variable the graph is with respect to Tij variations?
    key_dict, g = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path,
                                                                    edge_file_path=edge_file_path,
                                                                    key_dict_file_path=key_dict_file_path)
    if data_storage_location_keep_each_loop != None:
        if os.path.isfile(data_storage_location_keep_each_loop + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_num = last_row_explored["Iteration Number"].iloc[0]
        else:
            current_num = 0
            dictionary_fieldnames = ["Graph key", "Iteration Number", "objective_value", "total_capacity"]
            with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_num = 0

    objective_values = {}
    for num in range(current_num, num_iterations):
        for key in g.keys():
            new_Tij = {}
            for i, j in g[key].edges:
                if i < j:
                    new_Tij[(i, j)] =  (1 + np.random.normal(0, Tij_var))
            try:
                prob = cplex.Cplex()
                optimisation = Trusted_Node_Optimisation_No_Switching(prob, g[key], key_dict[key])
                sol_dict, prob, time_taken = optimisation.trusted_node_optimisation_program(time_limit=1e3,
                                                                                cdevices=c_source + c_detector,
                                                                                c_dist=c_fibre_per_km,
                                                                                c_conn_dist=c_initial_fibre_conn_per_km,
                                                                                cnetwork=overall_cost,
                                                                                Lambda=Lambda, Tij=new_Tij)
                flow_dict, binary_dict, lambda_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                objective_value = prob.solution.get_objective_value()
                if key == 2:
                    plot_graph(g[key], delta_dict=delta_dict, save_extra=f"Tij_iteration_{num}")
                total_capacity = 0.0
                for i,j in new_Tij.keys():
                    total_capacity += objective_value * new_Tij[(i,j)]
                if key not in objective_values.keys():
                    objective_values[key] = [(objective_value, total_capacity)]
                else:
                    objective_values[key].append((objective_value, total_capacity))

                if data_storage_location_keep_each_loop != None:
                    dictionary = [
                        {"Graph key": key, "Iteration Number": num,
                         "objective_value": objective_value, "total_capacity": total_capacity}]
                    dictionary_fieldnames = ["Graph key", "Iteration Number", "objective_value", "total_capacity"]
                    if os.path.isfile(data_storage_location_keep_each_loop + '.csv'):
                        with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                            writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                            writer.writerows(dictionary)
                    else:
                        with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                            writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                            writer.writeheader()
                            writer.writerows(dictionary)
            except:
                continue
    if data_storage_location_keep_each_loop != None:
        plot_information = pd.read_csv(data_storage_location_keep_each_loop + ".csv")
        for index, row in plot_information.iterrows():
            if row["Graph key"] not in objective_values.keys():
                objective_values[row["Graph key"]] = [(row["objective_value"], row["total_capacity"])]
            else:
                objective_values[row["Graph key"]].append((row["objective_value"], row["total_capacity"]))
    ####### WHAT TO DO HERE????
    ####### look at variation of objective and total capacity

    for key in objective_values.keys():
        objective_for_key = []
        total_cap_for_key = []
        objectives = objective_values[key]
        for i in range(len(objectives)):
            objective_for_key.append(objectives[i][0])
            total_cap_for_key.append(objectives[i][1])
        mean_objective = np.mean(objective_for_key)
        mean_total_cap = np.mean(total_cap_for_key)
        stdev_objective = np.std(objective_for_key)
        stdev_total_cap = np.std(total_cap_for_key)
        print("Mean for the objective is " + str(mean_objective) + ". StDev of objective is " + str(stdev_objective))
        print("Mean for the total capacity is " + str(mean_total_cap) + ". StDev of total capacity is " + str(stdev_total_cap))



def compare_no_switching_with_switching(c_source, c_detector, c_fibre_per_km, c_initial_fibre_conn_per_km, overall_cost_min, overall_cost_max,fswitch, Lambda, Tij, node_file_path, edge_file_path, key_dict_file_path,node_file_path_switching, edge_file_path_switching, key_dict_file_path_switching, data_storage_location_keep_each_loop = None, data_storage_location_keep_each_loop_switching = None):
    key_dict, g = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path,
                                                                    edge_file_path=edge_file_path,
                                                                    key_dict_file_path=key_dict_file_path)
    if data_storage_location_keep_each_loop != None:
        if os.path.isfile(data_storage_location_keep_each_loop + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_key = last_row_explored["Graph key"].iloc[0]
            current_cost = last_row_explored["overall_cost"].iloc[0]
        else:
            current_key = None
            current_cost = None
            dictionary_fieldnames = ["Graph key", "overall_cost", "objective_value"]
            with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_key = None
        current_cost = None

    objective_values = {}
    no_solution_list = []
    for overall_cost in np.arange(overall_cost_min, overall_cost_max, 10):
        if current_cost != None and current_cost != overall_cost:
            continue
        elif current_cost != None:
            current_cost = None
        for key in g.keys():
            if current_key != None and key != current_key:
                continue
            elif current_key != None:
                current_key = None
                continue
            if key not in no_solution_list:
                try:
                    prob = cplex.Cplex()
                    optimisation = Trusted_Node_Optimisation_No_Switching(prob, g[key], key_dict[key])
                    sol_dict, prob, time_taken = optimisation.trusted_node_optimisation_program(time_limit=2e2,
                                                                                    cdevices=c_source + c_detector,
                                                                                    c_dist=c_fibre_per_km,
                                                                                    c_conn_dist=c_initial_fibre_conn_per_km,
                                                                                    cnetwork=overall_cost,
                                                                                    Lambda=Lambda, Nmax = Lambda, Tij=Tij)
                    flow_dict, binary_dict, lambda_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                    objective_value = prob.solution.get_objective_value()
                    # if key == 2:
                    #     plot_graph(g[key], delta_dict=delta_dict, save_extra=f"comparison_no_switching_{key}_cost_{overall_cost}")
                    if key not in objective_values.keys():
                        objective_values[key] = [objective_value]
                    else:
                        objective_values[key].append(objective_value)

                    if data_storage_location_keep_each_loop != None:
                        dictionary = [
                            {"Graph key": key, "overall_cost": overall_cost,
                             "objective_value": objective_value}]
                        dictionary_fieldnames = ["Graph key","overall_cost", "objective_value"]
                        if os.path.isfile(data_storage_location_keep_each_loop + '.csv'):
                            with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                writer.writerows(dictionary)
                        else:
                            with open(data_storage_location_keep_each_loop + '.csv', mode='a') as csv_file:
                                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                writer.writeheader()
                                writer.writerows(dictionary)
                except:
                    continue
    key_dict_switching, g_switching = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path_switching,
                                                                    edge_file_path=edge_file_path_switching,
                                                                    key_dict_file_path=key_dict_file_path_switching)
    if data_storage_location_keep_each_loop_switching != None:
        if os.path.isfile(data_storage_location_keep_each_loop_switching + '.csv'):
            plot_information = pd.read_csv(data_storage_location_keep_each_loop_switching + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_key = last_row_explored["Graph key"].iloc[0]
            current_cost = last_row_explored["overall_cost"].iloc[0]
        else:
            current_key = None
            current_cost = None
            dictionary_fieldnames = ["Graph key","overall_cost", "objective_value"]
            with open(data_storage_location_keep_each_loop_switching + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_key = None
        current_cost = None

    objective_values_switching = {}
    no_solution_list_switching = []
    for overall_cost in np.arange(overall_cost_min, overall_cost_max, 10):
        if current_cost != None and current_cost != overall_cost:
            continue
        elif current_cost != None:
            current_cost = None
        for key in g_switching.keys():
            if current_key != None and key != current_key:
                continue
            elif current_key != None:
                current_key = None
                continue
            if key not in no_solution_list_switching:
                try:
                    key_dict_new = switched_fully_trusted_network.make_key_dict_bidirectional(key_dict[key])
                    prob = cplex.Cplex()
                    optimisation = switched_fully_trusted_network.Switched_Trusted_Node_Optimisation_Constant_Switching_Ratio(prob, g[key], key_dict_new)
                    sol_dict, prob = optimisation.trusted_node_optimisation_program(time_limit=2e2, csource=c_source, cdet=c_detector,
                                                                    c_dist=c_fibre_per_km, c_conn_dist=c_initial_fibre_conn_per_km, Lambda=Lambda, cgraph=overall_cost,
                                                                    f_switch=fswitch, Nmax = Lambda, Tij = Tij)
                    flow_dict, numbers_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
                    objective_value = prob.solution.get_objective_value()
                    # if key == 2:
                    #     plot_graph(g[key], delta_dict=delta_dict, save_extra=f"comparison_switching_{key}_cost_{overall_cost}")
                    if key not in objective_values_switching.keys():
                        objective_values_switching[key] = [objective_value]
                    else:
                        objective_values_switching[key].append(objective_value)

                    if data_storage_location_keep_each_loop_switching != None:
                        dictionary = [
                            {"Graph key": key, "overall_cost": overall_cost,
                             "objective_value": objective_value}]
                        dictionary_fieldnames = ["Graph key","overall_cost", "objective_value"]
                        if os.path.isfile(data_storage_location_keep_each_loop_switching + '.csv'):
                            with open(data_storage_location_keep_each_loop_switching + '.csv', mode='a') as csv_file:
                                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                writer.writerows(dictionary)
                        else:
                            with open(data_storage_location_keep_each_loop_switching + '.csv', mode='a') as csv_file:
                                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                                writer.writeheader()
                                writer.writerows(dictionary)
                except:
                    continue
    if data_storage_location_keep_each_loop != None:
        plot_information = pd.read_csv(data_storage_location_keep_each_loop + ".csv")
        for index, row in plot_information.iterrows():
            if row["overall_cost"] not in objective_values.keys():
                objective_values[row["overall_cost"]] = {row["Graph key"]: row["objective_value"]}
            else:
                objective_values[row["overall_cost"]][row["Graph key"]] = row["objective_value"]
    if data_storage_location_keep_each_loop_switching != None:
        plot_information = pd.read_csv(data_storage_location_keep_each_loop_switching + ".csv")
        for index, row in plot_information.iterrows():
            if row["overall_cost"] not in objective_values_switching.keys():
                objective_values_switching[row["overall_cost"]] = {row["Graph key"]: row["objective_value"]}
            else:
                objective_values_switching[row["overall_cost"]][row["Graph key"]] = row["objective_value"]

    objective_differences = {}
    for cost in objective_values.keys():
        for key in objective_values[cost].keys():
            if objective_values[cost][key] > 0.0001:
                if cost not in objective_differences.keys():
                    objective_differences[cost] = [(objective_values_switching[cost][key] - objective_values[cost][key]) / objective_values[cost][key]]
                else:
                    objective_differences[cost].append((objective_values_switching[cost][key] - objective_values[cost][key])/objective_values[cost][key])
    mean_objectives = {}
    std_objectives = {}
    for key in objective_differences.keys():
        mean_objectives[key] = np.mean(objective_differences[key])
        std_objectives[key] = np.std(objective_differences[key])
    mean_differences = []
    std_differences = []
    # topologies
    x = []
    for key in mean_objectives.keys():
        mean_differences.append(mean_objectives[key])
        std_differences.append(std_objectives[key])
        x.append(key)
    plt.errorbar(x, mean_differences, yerr=std_differences, color="r")
    plt.xlabel("Overall Cost of the Graph", fontsize=10)
    plt.ylabel("Percentage Difference of z_min between Switching and No Switching solutions", fontsize=10)
    # plt.legend(loc='upper right', fontsize='medium')
    plt.savefig("switching_no_switching_comparison.png")
    plt.show()

    objective_no_switching = {}
    for cost in objective_values.keys():
        for key in objective_values[cost].keys():
            if objective_values[overall_cost_min][key] > 0.0001:
                if cost not in objective_no_switching.keys():
                    objective_no_switching[cost] = [objective_values[cost][key] / objective_values[overall_cost_min][key]]
                else:
                    objective_no_switching[cost].append(objective_values[cost][key] /objective_values[overall_cost_min][key])
    mean_objectives = {}
    std_objectives = {}
    for key in objective_no_switching.keys():
        mean_objectives[key] = np.mean(objective_no_switching[key])
        std_objectives[key] = np.std(objective_no_switching[key])
    mean_differences = []
    std_differences = []
    # topologies
    x = []
    for key in mean_objectives.keys():
        mean_differences.append(mean_objectives[key])
        std_differences.append(std_objectives[key])
        x.append(key)
    plt.errorbar(x, mean_differences, yerr=std_differences, color="r")
    plt.xlabel("Overall Cost of the Graph", fontsize=10)
    plt.ylabel("z_min for current cost/z_min for minimum cost solution", fontsize=10)
    # plt.legend(loc='upper right', fontsize='medium')
    plt.savefig("no_switching_overall_cost_comparison.png")
    plt.show()


    objective_switching = {}
    for cost in objective_values_switching.keys():
        for key in objective_values_switching[cost].keys():
            if objective_values_switching[overall_cost_min][key] > 0.0001:
                if cost not in objective_switching.keys():
                    objective_switching[cost] = [objective_values_switching[cost][key] / objective_values_switching[overall_cost_min][key]]
                else:
                    objective_switching[cost].append(objective_values_switching[cost][key] /objective_values_switching[overall_cost_min][key])
    mean_objectives = {}
    std_objectives = {}
    for key in objective_switching.keys():
        mean_objectives[key] = np.mean(objective_switching[key])
        std_objectives[key] = np.std(objective_switching[key])
    mean_differences = []
    std_differences = []
    # topologies
    x = []
    for key in mean_objectives.keys():
        mean_differences.append(mean_objectives[key])
        std_differences.append(std_objectives[key])
        x.append(key)
    plt.errorbar(x, mean_differences, yerr=std_differences, color="r")
    plt.xlabel("Overall Cost of the Graph", fontsize=10)
    plt.ylabel("z_min for current cost/z_min for minimum cost solution", fontsize=10)
    # plt.legend(loc='upper right', fontsize='medium')
    plt.savefig("switching_overall_cost_comparison.png")
    plt.show()

#### Time taken investigations:

def time_variation_analysis(c_source, c_detector, c_fibre_per_km, c_initial_fibre_conn_per_km, overall_cost, Lambda, Tij, node_file_path, edge_file_path, key_dict_file_path, save_file):
    key_dict, g = utils.import_problem_trusted_nodes_given_key_dict(node_file_path=node_file_path,
                                                                    edge_file_path=edge_file_path,
                                                                    key_dict_file_path=key_dict_file_path)
    if save_file != None:
        if os.path.isfile(save_file + '.csv'):
            plot_information = pd.read_csv(save_file + ".csv")
            last_row_explored = plot_information.iloc[[-1]]
            current_key = last_row_explored["Graph key"].iloc[0]
            current_cost = last_row_explored["time_taken"].iloc[0]
        else:
            current_key = None
            current_cost = None
            dictionary_fieldnames = ["Graph key", "time_taken", "graph_size"]
            with open(save_file + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
    else:
        current_key = None
        current_cost = None

    time_taken = {}
    for key in g.keys():
        if current_key != None and current_key != key:
            continue
        elif current_key != None:
            current_key = None
            continue
        try:
            # key_dict_temp = trusted_nodes_utils.make_key_dict_bidirectional(key_dict[key])
            prob = cplex.Cplex()
            optimisation = Trusted_Node_Optimisation_No_Switching(prob, g[key], key_dict[key])
            sol_dict, prob, time = optimisation.trusted_node_optimisation_program(time_limit=3e3,
                                                                                        cdevices=c_source + c_detector,
                                                                                        c_dist=c_fibre_per_km,
                                                                                        c_conn_dist=c_initial_fibre_conn_per_km,
                                                                                        cnetwork=overall_cost,
                                                                                        Lambda=Lambda, Nmax=Lambda,
                                                                                        Tij=Tij)
            # optim = Optimisation_Switching_Calibration_fixed_frac_calibration_time(prob=prob, g=g[key],
            #                                                                        key_dict=key_dict_temp)
            # sol_dict, prob, time_taken = optim.initial_optimisation_cost_reduction(cmin = cmin, time_limit=time_limit, cost_on_trusted_node=cost_on_trusted_node,
            #                                            cost_detector=cost_detector, cost_source=cost_source, f_switch=f_switch, Lambda=Lambda)
            number_nodes = g[key].number_of_nodes()
            if number_nodes in time_taken.keys():
                time_taken[number_nodes].append(time)
            else:
                time_taken[number_nodes] = [time]
            if save_file != None:
                dictionary = [
                    {"Graph key": key,"time_taken": time, "graph_size" : number_nodes}]
                dictionary_fieldnames = ["Graph key", "time_taken", "graph_size"]
                if os.path.isfile(save_file + '.csv'):
                    with open(save_file + '.csv', mode='a') as csv_file:
                        writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                        writer.writerows(dictionary)
                else:
                    with open(save_file + '.csv', mode='a') as csv_file:
                        writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                        writer.writeheader()
                        writer.writerows(dictionary)
        except:
            continue

    if save_file != None:
        plot_information = pd.read_csv(save_file + ".csv")
        for index, row in plot_information.iterrows():
            if row["graph_size"] not in time_taken.keys():
                time_taken[row["graph_size"]] = [row["time_taken"]]
            else:
                time_taken[row["graph_size"]].append(row["time_taken"])
    time_costs_mean_std = {}
    x = []
    y = []
    yerr = []
    for key in time_taken:
        if key < 14:

            time_costs_mean_std[key] = [np.mean(time_taken[key]), np.std(time_taken[key])]
            x.append(key)
            y.append(time_costs_mean_std[key][0])
            yerr.append(time_costs_mean_std[key][1])
        # x                     y               yerr
    # fit of exponential curve to initial points
    # x_exponential = x
    # y_exponential = y
    # popt, pcov = curve_fit(optimisation_switched.exponential_fit, x_exponential, y_exponential)
    # x_exponential = np.arange(x[0], x[-1], 0.1)
    # y_exponential = [optimisation_switched.exponential_fit(a, popt[0], popt[1], popt[2]) for a in x_exponential]
    # # fit of polynomial
    # popt_poly, pcov_poly = curve_fit(optimisation_switched.polynomial, x, y)
    # y_poly = [optimisation_switched.polynomial(a, popt_poly[0], popt_poly[1]) for a in x_exponential]

    plt.errorbar(x, y, yerr=yerr, color="r")
    plt.yscale(value = "log")
    # plt.plot(x_exponential[:int(np.ceil(len(x_exponential)/1.25))], y_exponential[:int(np.ceil(len(x_exponential)/1.25))], color = "b")
    # plt.plot(x_exponential, y_poly, color = "k")
    plt.xlabel("Number of Nodes in Graph", fontsize=10)
    plt.ylabel("Time/s", fontsize=10)
    # plt.legend(loc='upper right', fontsize='medium')
    plt.savefig("time_investigation_mesh_topology.png")
    plt.show()


### Investigation of Optimal Graphs

# def plot_optimal_graphs_for_different_costs():

 #### Standard price for source: 1.5k, Standard price for detectors 120k, per fibre km: 3.8k, init cost to bury: 32.7k
#### if we set the price of source to 0.01 --> det price: 0.8, per fibre km: 0.025, init cost to bury: 0.22 (set this to 0.2)

if __name__ == "__main__":
    # compare_no_switching_with_switching(c_source = 0.01, c_detector= 0.8, c_fibre_per_km = 0.025, c_initial_fibre_conn_per_km = 0.22,
    #                                         overall_cost_min = 100, overall_cost_max= 400, fswitch= 0.1, Lambda = 100, Tij = None, node_file_path = "1_trusted_nodes_positions.csv", edge_file_path = "1_trusted_nodes_capacities.csv", key_dict_file_path = "1_trusted_nodes_key_dict.csv",
    #                                         node_file_path_switching = "1_trusted_nodes_positions.csv",
    #                                         edge_file_path_switching = "1_trusted_nodes_capacities.csv", key_dict_file_path_switching = "1_trusted_nodes_key_dict.csv")

    cost_detector_investigation_switching(c_source = 0.01, c_detector = 1.5, c_fibre_per_km = 0.025, c_initial_fibre_conn_per_km = 0,
                                          overall_cost = 20, frac_switch = 0.1, Lambda = 100, Tij = None, node_file_path = "3_trusted_nodes_positions.csv", edge_file_path = "3_trusted_nodes_capacities_switched.csv",
                                          edge_no_switching_path= "3_trusted_nodes_capacities.csv",
                                          key_dict_file_path = "3_trusted_nodes_key_dict.csv", data_storage_location_keep_each_loop_no_switching = "switching_detector_cost_variation_no_switching_0_init_fibre_cost_20_total_cost_new",
                                          data_storage_location_keep_each_loop_switching ="switching_detector_cost_variation_switching_0_init_fibre_cost_20_total_cost_new")


    # cost_overall_investigation_switching(c_source = 0.01, c_detector = 0.8, c_fibre_per_km = 0.025, c_initial_fibre_conn_per_km = 0.25,
    #                                       overall_cost = 300, frac_switch = 0.1, Lambda = 100, Tij = None, node_file_path = "3_trusted_nodes_positions.csv", edge_file_path = "3_trusted_nodes_capacities_switched.csv",
    #                                       edge_no_switching_path= "3_trusted_nodes_capacities.csv",
    #                                       key_dict_file_path = "3_trusted_nodes_key_dict.csv", data_storage_location_keep_each_loop_no_switching = "switching_overall_cost_variation_no_switching_corr",
    #                                       data_storage_location_keep_each_loop_switching ="switching_overall_cost_variation_switching_corr")


    # cost_investigations(c_source = 0.01, c_detector= 1, c_fibre_per_km = 0.25, c_initial_fibre_conn_per_km= 1.0, overall_cost = 400, Lambda = 100, Tij = None,
    #                     node_file_path = "2_trusted_nodes_positions.csv", edge_file_path = "2_trusted_nodes_capacities.csv", key_dict_file_path = "2_trusted_nodes_key_dict.csv", data_storage_location_keep_each_loop = "cost_analysis_no_switching")
    # cost_investigations_shorter(c_source = 0.01, c_detector= 1, c_fibre_per_km = 0.2, c_initial_fibre_conn_per_km= 1.0, overall_cost = 400, Lambda = 100, Tij = None,
    #                     node_file_path = "2_trusted_nodes_positions.csv", edge_file_path = "2_trusted_nodes_capacities.csv", key_dict_file_path = "2_trusted_nodes_key_dict.csv",
    #                             data_storage_location_keep_each_loop_fibre_init_fibre = "cost_fibre_init_fibre.csv",
    #                             data_storage_location_keep_each_loop_detect_init_fibre = "cost_detector_init_fibre.csv",
    #                             data_storage_location_keep_each_loop_detect_fibre = "cost_detector_fibre.csv")
    # variation_of_fibre_cost_effects_with_number_channels(c_source = 0.01, c_detector =0.8, c_fibre_per_km = 0.2,
    #                                                      c_initial_fibre_conn_per_km = 1.0, overall_cost = 400, Lambda = 100, Tij = None,
    #                                                      node_file_path="2_trusted_nodes_positions.csv",
    #                                                      edge_file_path="2_trusted_nodes_capacities.csv",
    #                                                      key_dict_file_path="2_trusted_nodes_key_dict.csv",
    #                                                      data_storage_location_keep_each_loop_fibre_init_fibre="cost_fibre_different_no_connections_corrected.csv",
    #                                                      data_storage_location_keep_each_loop_init_fibre = "file_to_delete.csv")

    # compare_no_switching_with_switching(c_source = 0.02, c_detector = 0.8, c_fibre_per_km = 0.025, c_initial_fibre_conn_per_km = 0.2,
    #                                     overall_cost_min = 300, overall_cost_max = 510, fswitch = 0.1, Lambda = 100, Tij = None, node_file_path="3_trusted_nodes_positions.csv",
    #                                                      edge_file_path="3_trusted_nodes_capacities.csv",
    #                                                      key_dict_file_path="3_trusted_nodes_key_dict.csv", node_file_path_switching = "3_trusted_nodes_positions.csv",
    #                                     edge_file_path_switching = "3_trusted_nodes_capacities_switched.csv", key_dict_file_path_switching = "3_trusted_nodes_key_dict.csv",
    #                                     data_storage_location_keep_each_loop="switching_no_switching_comparison_no_switching_data.csv",
    #                                     data_storage_location_keep_each_loop_switching="switching_no_switching_comparison_switching_data.csv")

    # comparison_switch_no_switch_different_no_channels(c_source = 0.02, c_detector = 0.8, c_fibre_per_km = 0.025, c_initial_fibre_conn_per_km = 0.2,
    #                                     overall_cost_min = 300, overall_cost_max = 510, fswitch = 0.1, Lambda = 100, Tij = None, node_file_path="3_trusted_nodes_positions.csv",
    #                                                      edge_file_path="3_trusted_nodes_capacities.csv",
    #                                                      key_dict_file_path="3_trusted_nodes_key_dict.csv", node_file_path_switching = "3_trusted_nodes_positions.csv",
    #                                     edge_file_path_switching = "3_trusted_nodes_capacities_switched.csv", key_dict_file_path_switching = "3_trusted_nodes_key_dict.csv",
    #                                     data_storage_location_keep_each_loop="switching_no_switching_comparison_different_fibre_channels_no_switching_data.csv",
    #                                     data_storage_location_keep_each_loop_switching="switching_no_switching_comparison_different_fibre_channels_switching_data_corr.csv")

    # l_switch_parameter_sweep(node_file_path = "4_trusted_nodes_positions.csv", edge_file_path = "4_trusted_nodes_capacities.csv", key_dict_file_path = "4_trusted_nodes_key_dict.csv", node_file_path_switching = "4_trusted_nodes_positions.csv",
    #                          edge_file_path_switching = "4_trusted_nodes_capacities_switched_db_", key_dict_file_path_switching = "4_trusted_nodes_key_dict.csv", overall_cost = 300,
    #                          c_source = 0.02, f_switch = 0.1, Tij=None, c_detector=0.8, c_fibre_per_km=0.025,
    #                          c_initial_fibre_conn_per_km=0.2, Lambda=100, data_storage_location_keep_each_loop="switching_loss_data_switching",
    #                          data_storage_location_keep_each_loop_no_switch="switching_loss_data_no_switching")

    # f_switch_parameter_sweep(c_source=0.02, c_detector=0.8, c_fibre_per_km=0.025,
    #                                     c_initial_fibre_conn_per_km=0.2,
    #                                     overall_cost=300, Lambda=100, Tij=None,
    #                                     node_file_path="3_trusted_nodes_positions.csv",
    #                                     edge_file_path="3_trusted_nodes_capacities.csv",
    #                                     key_dict_file_path="3_trusted_nodes_key_dict.csv",
    #                                     node_file_path_switching="3_trusted_nodes_positions.csv",
    #                                     edge_file_path_switching="3_trusted_nodes_capacities_switched.csv",
    #                                     key_dict_file_path_switching="3_trusted_nodes_key_dict.csv",
    #                                     data_storage_location_keep_each_loop = "switching_no_switching_fswitch_switching_corrected", data_storage_location_keep_each_loop_no_switch = "switching_no_switching_fswitch_no_switching_corrected")

    # time_variation_analysis(c_source= 0.01, c_detector = 0.8, c_fibre_per_km = 0.025, c_initial_fibre_conn_per_km = 0.22, overall_cost = 400, Lambda =100,
    #                         Tij = None, node_file_path = "1_trusted_nodes_positions.csv", edge_file_path = "1_trusted_nodes_capacities.csv", key_dict_file_path = "1_trusted_nodes_key_dict.csv", save_file = "time_taken_investigation")
    # key_dict, g = utils.import_problem_trusted_nodes_given_key_dict(node_file_path="small_node_graph_test_nodes_1.csv",edge_file_path = "small_node_graph_test_1.csv", key_dict_file_path="small_node_key_dict.csv")
    # for key in g.keys():
    #     prob = cplex.Cplex()
    #     optimisation = Trusted_Node_Optimisation_No_Switching(prob, g[key], key_dict[key])
    #     sol_dict, prob, time_taken = optimisation.trusted_node_optimisation_program(Lambda =10,)
    #     flow_dict, binary_dict, lambda_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
    #     plot_graph(g[key], delta_dict)
    #     for i in sol_dict.keys():
    #         print(str(i) + ": " + str(sol_dict[i]))




