import cplex
import networkx as nx
import csv
import utils
import time
import matplotlib.pyplot as plt


class Switched_Trusted_Node_Optimisation():

    def __init__(self, prob, g, key_dict):
        self.prob = prob
        self.g = g
        self.key_dict = key_dict
        super().__init__()

    def log_optimal_solution_to_problem(self, save_file, graph_id):
        sol_dict = self.create_sol_dict(self.prob)
        use_dict, number_dict, delta_dict = self.split_sol_dict(sol_dict)
        ## if file does not exist - we wish to store information of q_{i,j,d}^{m}, w_{i,j,d}^{m}, lambda_{d}^{m}, delta_{i,j,d}^{m}
        ## need to generate all possible values of i,j,d available. Need to think of the most appropriate way to store this data.
        dict = {"ID": graph_id}
        dict.update(use_dict)
        dict.update(number_dict)
        dict.update(delta_dict)
        dictionary = [dict]
        dictionary_fieldnames = list(dict.keys())
        with open(save_file + '.csv', mode='a') as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
            writer.writeheader()
            writer.writerows(dictionary)



    def split_sol_dict(self, sol_dict):
        """
        Split the solution dictionary into 2 dictionaries containing the fractional usage variables only and the binary
        variables only
        Parameters
        ----------
        sol_dict : The solution dictionary containing solutions to the primary flow problem

        Returns : A dictionary with only the fractional detectors used, and a dictionary with only the binary values of
                whether the detector is on or off for cold and hot
        -------

        """
        use_dict = {}
        number_dict = {}
        delta_dict = {}
        for key in sol_dict:
            # get all keys that are flow and add to dictionary
            if key[0] == "x" or key[0] == "z":
                use_dict[key] = sol_dict[key]
            elif key[0] == "N":
                number_dict[key] = sol_dict[key]
            elif key[0] == "d":
                delta_dict[key] = sol_dict[key]
        return use_dict, number_dict, delta_dict


    def create_sol_dict(self):
        """
        Create a dictionary with the solution of the parameters
        """
        names = self.prob.variables.get_names()
        values = self.prob.solution.get_values()
        sol_dict = {names[idx]: (values[idx]) for idx in range(self.prob.variables.get_num())}
        return sol_dict

    def flow_conservation_constraint(self):
        pass

    def no_flow_into_source_out_sink(self):
        pass

    def cost_constraint(self, *args, **kwargs):
        pass

    def add_source_capacity_on_edge_constraint(self):
        pass

    def add_detector_capacity_on_edge_constraint(self):
        pass

    def fibre_constraint(self):
        pass

    def add_fibre_on_constraint(self):
        pass

    def maximise_lowest_capacity_connection_objective(self, Tij):
        pass



class Switched_Trusted_Node_Optimisation_Switching_W_terms(Switched_Trusted_Node_Optimisation):

    def __init__(self, prob, g, key_dict):
        super().__init__(prob = prob, g = g, key_dict = key_dict)

    def flow_conservation_constraint(self):
        """
        Constraint for flow conservation: \sum_{m \in N(i)} x^{k}_{(i,m)} + x^{k_R}_{(m,i)} - x^{k}_{(m,i)} - x^{k_R}_{(i,m)} = 0
        """
        variable_names = [f'x{i}_{j}_k{k[0]}_{k[1]}' for k in self.key_dict for i, j in list(self.g.edges)]
        self.prob.variables.add(names=variable_names, types=[self.prob.variables.type.continuous] * len(variable_names))
        for i in self.g.nodes:
            for k in self.key_dict:
                if k[0] < k[1] and k[0] != i and k[1] != i:
                    variables = []
                    val = []
                    for n in self.g.neighbors(i):
                        variables.extend(
                            [f"x{i}_{n}_k{k[0]}_{k[1]}", f"x{n}_{i}_k{k[1]}_{k[0]}", f"x{n}_{i}_k{k[0]}_{k[1]}",
                             f"x{i}_{n}_k{k[1]}_{k[0]}"])
                        val.extend([1, 1, -1, -1])
                    lin_expressions = [cplex.SparsePair(ind=variables, val=val)]
                    self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["E"], rhs=[0.])

    def no_flow_into_source_out_sink(self):
        """
        Constraint for the prevention of flow into the source and out of the sink:
                x^{k=(i,j)}_{(j,m)} + x^{k_R=(j,i)}_{(m,j)} = 0
                x^{k=(i,j)}_{(m,i)} + x^{k_R=(j,i)}_{(i,m)} = 0
        """
        for k in self.key_dict:
            if k[0] < k[1]:
                for n in self.g.neighbors(k[1]):
                    variables = [f"x{k[1]}_{n}_k{k[0]}_{k[1]}", f"x{n}_{k[1]}_k{k[1]}_{k[0]}"]
                    coeffs = [1, 1]
                    lin_expressions = [cplex.SparsePair(ind=variables, val=coeffs)]
                    self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["E"], rhs=[0.])
                for n in self.g.neighbors(k[0]):
                    variables = [f"x{n}_{k[0]}_k{k[0]}_{k[1]}", f"x{k[0]}_{n}_k{k[1]}_{k[0]}"]
                    coeffs = [1, 1]
                    lin_expressions = [cplex.SparsePair(ind=variables, val=coeffs)]
                    self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["E"], rhs=[0.])

    def cost_constraint(self, csource, cdet, cfibres, cconn, cgraph, max_devices, * args, **kwargs):
        """
        Adds the constraint to restrict the cost of the network:
        \sum_{i \in V} C_{source}N_{i}^{S}+ C_{det}N_{i}^{D} + \sum_{(i,j) \in E} C_{i,j}^{fibre}N_{i,j}^{fibre} + C_{i,j}^{conn}
        \delta_{i,j} \leq C

        """
        # cfibre is a dictionary with fibre costs based on the edge
        # NEED A FUNCTION TO GET CFIBRE
        source_terms = []
        detector_terms = []
        for n in self.g.nodes:
            source_terms.append(f"N_{n}_S")
            detector_terms.append(f"N_{n}_D")
        self.prob.variables.add(names=source_terms, types=[self.prob.variables.type.integer] * len(source_terms),
                           ub=[max_devices] * len(source_terms))
        self.prob.variables.add(names=detector_terms, types=[self.prob.variables.type.integer] * len(detector_terms),
                           ub=[max_devices] * len(detector_terms))
        delta_terms = []
        fibre_terms = []
        for i, j in self.g.edges:
            if i < j:
                fibre_terms.append(f"N_{i}_{j}_fibre")
                delta_terms.append(f"delta_{i}_{j}")
        self.prob.variables.add(names=fibre_terms, types=[self.prob.variables.type.integer] * len(fibre_terms))
        self.prob.variables.add(names=delta_terms, types=[self.prob.variables.type.binary] * len(delta_terms))
        variables = []
        val = []
        for i in self.g.nodes:
            variables.extend([f"N_{i}_S", f"N_{i}_D"])
            val.extend([csource, cdet])
        for i, j in self.g.edges:
            if i < j:
                variables.extend([f"N_{i}_{j}_fibre", f"delta_{i}_{j}"])
                val.extend([cfibres[(i, j)], cconn[(i, j)]])

        lin_expressions = [cplex.SparsePair(ind=variables, val=val)]
        self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["L"], rhs=[cgraph])

    def add_untransformed_transformed_relationship(self, NT):
        """
        Adds the constraint that relates the time after considering switching time with the time before switching:
        \sum_{k \in K} x_{i,j}^{k} = \sum_{k \in K} w_{i,j}^{k} - \frac{N}{T}N_{j}^{D}
        """
        variable_names = [f'w{i}_{j}_k{k[0]}_{k[1]}' for k in self.key_dict for i, j in list(self.g.edges)]
        self.prob.variables.add(names=variable_names, types=[self.prob.variables.type.continuous] * len(variable_names))
        for i, j in self.g.edges:
            ind = []
            val = []
            for k in self.key_dict:
                ind.append(f"x{i}_{j}_k{k[0]}_{k[1]}")
                val.append(1)
                ind.append(f"w{i}_{j}_k{k[0]}_{k[1]}")
                val.append(-1)
            ind.append(f"N_{j}_D")
            val.append(NT)
            lin_expressions = [cplex.SparsePair(ind=ind, val=val)]
            self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=['E'], rhs=[0.])

    def add_source_capacity_on_edge_constraint(self):
        """
        Adds the constraint to restrict the flow along an edge to the number of sources at the source node:
        \sum_{j \in N(i)} \frac{\sum_{k \in K} w_{i,j}^{k}}{c_{i,j}} \leq N^S_{i}
        """
        for i in self.g.nodes:
            flow = []
            coeff = []
            for j in self.g.neighbors(i):
                capacity = self.g.edges[[i, j]]["capacity"]
                if capacity > 0.000001:
                    for k in self.key_dict:
                        flow.append(f"w{i}_{j}_k{k[0]}_{k[1]}")
                        coeff.append(1 / capacity)
                else:
                    for k in self.key_dict:
                        self.prob.linear_constraints.add(
                            lin_expr=[cplex.SparsePair(ind=[f"w{i}_{j}_k{k[0]}_{k[1]}"], val=[1])], senses="E",
                            rhs=[0.])
            flow.append(f"N_{i}_S")
            coeff.append(-1)
            lin_expressions = [cplex.SparsePair(ind=flow, val=coeff)]
            self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=['L'], rhs=[0.])

    def add_detector_capacity_on_edge_constraint(self):
        """
        Adds the constraint to restrict the flow along an edge to the number of detectors at the sink:
        \sum_{j \in N(i)} \frac{\sum_{k \in K} w_{j,i}^{k}}{c_{i,j}} \leq N^D_{i}
        """
        for i in self.g.nodes:
            flow_along_edge = []
            coeff = []
            for j in self.g.adj[i]:
                capacity = int(self.g.edges[[i, j]]["capacity"])
                if capacity > 0.0000001:
                    for k in self.key_dict:
                        flow_along_edge.append(f"w{j}_{i}_k{k[0]}_{k[1]}")
                        coeff.append(1 / capacity)
                else:
                    for k in self.key_dict:
                        self.prob.linear_constraints.add(
                            lin_expr=[cplex.SparsePair(ind=[f"w{j}_{i}_k{k[0]}_{k[1]}"], val=[1])],
                            senses="E", rhs=[0.])

            flow_along_edge.append(f"N_{i}_D")
            coeff.append(-1)
            lin_expressions = [cplex.SparsePair(ind=flow_along_edge, val=coeff)]
            self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=['L'], rhs=[0.])

    def fibre_constraint(self):
        """
        Adds the fibre constraint to the problem - flow across an edge cannot exceed fibre capacity:
        \sum_{k \in K} w_{(i,j)}^{k} + w_{(j,i)}^{k} \leq c_{i,j}(N_{i,j}^{fibre})

        """
        for i, j in self.g.edges:
            if i < j:
                variables = []
                val = []
                for k in self.key_dict:
                    variables.extend([f"w{i}_{j}_k{k[0]}_{k[1]}", f"w{j}_{i}_k{k[0]}_{k[1]}"])
                    val.extend([1, 1])
                variables.append(f"N_{i}_{j}_fibre")
                val.append(-int(self.g.edges[[i, j]]["capacity"]))
                lin_expressions = [cplex.SparsePair(ind=variables, val=val)]
                self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["L"], rhs=[0.])

    def source_more_than_detectors_approx(self):
        """
        Adds the constraint that approximates the requirement for more sources than detectors on every adjacent node:
        N_{j}^{D} \leq N_{i}^{S} \forall j \in N(i)
        """
        for i in self.g.nodes:
            for j in self.g.adj[i]:
                variables = [f"N_{j}_D", f"N_{i}_S"]
                val = [1, -1]
                lin_expressions = [cplex.SparsePair(ind=variables, val=val)]
                self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["L"], rhs=[0.])

    def add_fibre_on_constraint(self, Lambda):
        """
        adds the constraint that a fibre must be on to have any fibres N_{i,j}^{fibres} \leq Lambda \delta_{i,j}
        """
        for i, j in self.g.edges:
            if i < j:
                variables = [f"N_{i}_{j}_fibre", f"delta_{i}_{j}"]
                val = [1, -Lambda]
                lin_expressions = [cplex.SparsePair(ind=variables, val=val)]
                self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["L"], rhs=[0.])

    def maximise_lowest_capacity_connection_objective(self, Tij):
        """
        Adds the maximise lowest capacity objective:
        maximise z
        subject to:
        z \leq \sum_{n \in N(j)} x_{(n,j)}^{k=(i,j)} + x_{(j,n)}^{k=(j,i)}
        """
        variable_names = ['z']
        self.prob.variables.add(names=variable_names, types=[self.prob.variables.type.continuous])
        for k in self.key_dict:
            variables = []
            val = []
            for n in self.g.neighbors(k[1]):

                if Tij == None:
                    variables.extend([f"x{n}_{k[1]}_k{k[0]}_{k[1]}", f"x{k[1]}_{n}_k{k[1]}_{k[0]}"])
                    val.append(-1)
                    val.append(-1)
                elif Tij[k] > 0.000001:
                    variables.extend([f"x{n}_{k[1]}_k{k[0]}_{k[1]}", f"x{k[1]}_{n}_k{k[1]}_{k[0]}"])
                    val.append(-1 / Tij[k])
                    val.append(-1 / Tij[k])
            variables.append("z")
            val.append(1)
            lin_expressions = [cplex.SparsePair(ind=variables, val=val)]
            self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["L"], rhs=[0.])
        obj_vals = [("z", 1)]
        self.prob.objective.set_linear(obj_vals)
        self.prob.objective.set_sense(self.prob.objective.sense.maximize)

    def trusted_node_optimisation_program(self, time_limit=1e5, csource=0.1, cdet=1, c_dist=0.01,
                                          c_conn_dist=0.1, Lambda=10, cgraph=100, Nmax=10, NT=0, Tij = None):
        """
        set up and solve the problem for minimising the overall cost of the network
        """
        print("Start Optimisation")
        self.flow_conservation_constraint()
        # add_capacity_constraint(prob, g, key_dict, Lambda=Lambda)
        self.no_flow_into_source_out_sink()
        cfibres, cconn = get_cfibre(self.g, c_dist, c_conn_dist)
        self.cost_constraint(csource, cdet, cfibres, cconn, cgraph, Nmax)
        # add_minimise_trusted_nodes_objective(prob, g)
        self.add_untransformed_transformed_relationship(NT)
        self.add_source_capacity_on_edge_constraint()
        self.add_detector_capacity_on_edge_constraint()
        self.fibre_constraint()
        self.source_more_than_detectors_approx()
        self.add_fibre_on_constraint(Lambda)
        self.maximise_lowest_capacity_connection_objective(Tij = Tij)

        self.prob.write("test_multi.lp")
        self.prob.parameters.lpmethod.set(3)
        self.prob.parameters.mip.limits.cutpasses.set(1)
        self.prob.parameters.mip.strategy.probe.set(-1)
        self.prob.parameters.mip.strategy.variableselect.set(2)
        self.prob.parameters.mip.strategy.kappastats.set(1)
        self.prob.parameters.mip.tolerances.mipgap.set(float(0.01))
        # prob.parameters.simplex.limits.iterations = 50
        print(self.prob.parameters.get_changed())
        self.prob.parameters.timelimit.set(time_limit)
        t_1 = time.time()

        self.prob.solve()
        t_2 = time.time()
        print("Time to solve problem: " + str(t_2 - t_1))
        print(f"The minimum Cost of Network: {self.prob.solution.get_objective_value()}")
        print(f"Number of Variables = {self.prob.variables.get_num()}")
        print(f"Number of Conditions = {self.prob.linear_constraints.get_num()}")
        sol_dict = self.create_sol_dict()
        flow_dict, numbers_dict, delta_dict = self.split_sol_dict(sol_dict)
        return sol_dict, self.prob


class Switched_Trusted_Node_Optimisation_Constant_Switching_Ratio(Switched_Trusted_Node_Optimisation):


    def __init__(self, prob, g, key_dict):
        super().__init__(prob = prob, g = g, key_dict = key_dict)

    def flow_conservation_constraint(self):
        """
        Constraint for flow conservation: \sum_{m \in N(i)} x^{k}_{(i,m)} + x^{k_R}_{(m,i)} - x^{k}_{(m,i)} - x^{k_R}_{(i,m)} = 0
        """
        variable_names = [f'x{i}_{j}_k{k[0]}_{k[1]}' for k in self.key_dict for i, j in list(self.g.edges)]
        self.prob.variables.add(names=variable_names, types=[self.prob.variables.type.continuous] * len(variable_names))
        for i in self.g.nodes:
            for k in self.key_dict:
                if k[0] < k[1] and k[0] != i and k[1] != i:
                    variables = []
                    val = []
                    for n in self.g.neighbors(i):
                        variables.extend(
                            [f"x{i}_{n}_k{k[0]}_{k[1]}", f"x{n}_{i}_k{k[1]}_{k[0]}", f"x{n}_{i}_k{k[0]}_{k[1]}",
                             f"x{i}_{n}_k{k[1]}_{k[0]}"])
                        val.extend([1, 1, -1, -1])
                    lin_expressions = [cplex.SparsePair(ind=variables, val=val)]
                    self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["E"], rhs=[0.])

    def no_flow_into_source_out_sink(self):
        """
        Constraint for the prevention of flow into the source and out of the sink:
                x^{k=(i,j)}_{(j,m)} + x^{k_R=(j,i)}_{(m,j)} = 0
                x^{k=(i,j)}_{(m,i)} + x^{k_R=(j,i)}_{(i,m)} = 0
        """
        for k in self.key_dict:
            if k[0] < k[1]:
                for n in self.g.neighbors(k[1]):
                    variables = [f"x{k[1]}_{n}_k{k[0]}_{k[1]}", f"x{n}_{k[1]}_k{k[1]}_{k[0]}"]
                    coeffs = [1, 1]
                    lin_expressions = [cplex.SparsePair(ind=variables, val=coeffs)]
                    self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["E"], rhs=[0.])
                for n in self.g.neighbors(k[0]):
                    variables = [f"x{n}_{k[0]}_k{k[0]}_{k[1]}", f"x{k[0]}_{n}_k{k[1]}_{k[0]}"]
                    coeffs = [1, 1]
                    lin_expressions = [cplex.SparsePair(ind=variables, val=coeffs)]
                    self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["E"], rhs=[0.])

    def cost_constraint(self, csource, cdet, cfibres, cconn, cgraph, max_devices, * args, **kwargs):
        """
        Adds the constraint to restrict the cost of the network:
        \sum_{i \in V} C_{source}N_{i}^{S}+ C_{det}N_{i}^{D} + \sum_{(i,j) \in E} C_{i,j}^{fibre}N_{i,j}^{fibre} + C_{i,j}^{conn}
        \delta_{i,j} \leq C

        """
        # cfibre is a dictionary with fibre costs based on the edge
        # NEED A FUNCTION TO GET CFIBRE
        source_terms = []
        detector_terms = []
        ub = []
        for n in self.g.nodes:
            ub.append(max_devices * len(list(self.g.neighbors(n))))
            source_terms.append(f"N_{n}_S")
            detector_terms.append(f"N_{n}_D")
        self.prob.variables.add(names=source_terms, types=[self.prob.variables.type.integer] * len(source_terms),
                           ub=ub)
        self.prob.variables.add(names=detector_terms, types=[self.prob.variables.type.integer] * len(detector_terms),
                           ub=ub)
        delta_terms = []
        fibre_terms = []
        for i, j in self.g.edges:
            if i < j:
                fibre_terms.append(f"N_{i}_{j}_fibre")
                delta_terms.append(f"delta_{i}_{j}")
        self.prob.variables.add(names=fibre_terms, types=[self.prob.variables.type.integer] * len(fibre_terms))
        self.prob.variables.add(names=delta_terms, types=[self.prob.variables.type.binary] * len(delta_terms))
        variables = []
        val = []
        for i in self.g.nodes:
            variables.extend([f"N_{i}_S", f"N_{i}_D"])
            val.extend([csource, cdet])
        for i, j in self.g.edges:
            if i < j:
                variables.extend([f"N_{i}_{j}_fibre", f"delta_{i}_{j}"])
                val.extend([cfibres[(i, j)], cconn[(i, j)]])

        lin_expressions = [cplex.SparsePair(ind=variables, val=val)]
        self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["L"], rhs=[float(cgraph)])

    def add_source_capacity_on_edge_constraint(self, f_switch):
        """
        Adds the constraint to restrict the flow along an edge to the number of sources at the source node:
        \sum_{j \in N(i)} \frac{\sum_{k \in K} x_{i,j}^{k}}{c_{i,j}} \leq (1-f_switch) N^S_{i}
        """
        for i in self.g.nodes:
            flow = []
            coeff = []
            for j in self.g.neighbors(i):
                capacity = self.g.edges[[i, j]]["capacity"]
                if capacity > 0.000001:
                    for k in self.key_dict:
                        flow.append(f"x{i}_{j}_k{k[0]}_{k[1]}")
                        coeff.append(1 / capacity)
                else:
                    for k in self.key_dict:
                        self.prob.linear_constraints.add(
                            lin_expr=[cplex.SparsePair(ind=[f"x{i}_{j}_k{k[0]}_{k[1]}"], val=[1])], senses="E",
                            rhs=[0.])
            flow.append(f"N_{i}_S")
            coeff.append(-(1 - f_switch))
            lin_expressions = [cplex.SparsePair(ind=flow, val=coeff)]
            self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=['L'], rhs=[0.])

    def add_detector_capacity_on_edge_constraint(self, f_switch):
        """
        Adds the constraint to restrict the flow along an edge to the number of detectors at the sink:
        \sum_{j \in N(i)} \frac{\sum_{k \in K} w_{j,i}^{k}}{c_{i,j}} \leq (1-f_switch)N^D_{i}
        """
        for i in self.g.nodes:
            flow_along_edge = []
            coeff = []
            for j in self.g.adj[i]:
                capacity = int(self.g.edges[[i, j]]["capacity"])
                if capacity > 0.0000001:
                    for k in self.key_dict:
                        flow_along_edge.append(f"x{j}_{i}_k{k[0]}_{k[1]}")
                        coeff.append(1 / capacity)
                else:
                    for k in self.key_dict:
                        self.prob.linear_constraints.add(
                            lin_expr=[cplex.SparsePair(ind=[f"x{j}_{i}_k{k[0]}_{k[1]}"], val=[1])],
                            senses="E", rhs=[0.])

            flow_along_edge.append(f"N_{i}_D")
            coeff.append(-(1 - f_switch))
            lin_expressions = [cplex.SparsePair(ind=flow_along_edge, val=coeff)]
            self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=['L'], rhs=[0.])

    def fibre_constraint(self, f_switch, no_channels = 1):
        """
        Adds the fibre constraint to the problem - flow across an edge cannot exceed fibre capacity:
        \sum_{k \in K} w_{(i,j)}^{k} + w_{(j,i)}^{k} \leq N_channs(1-f_switch)c_{i,j}(N_{i,j}^{fibre})

        """
        for i, j in self.g.edges:
            if i < j:
                variables = []
                val = []
                for k in self.key_dict:
                    variables.extend([f"x{i}_{j}_k{k[0]}_{k[1]}", f"x{j}_{i}_k{k[0]}_{k[1]}"])
                    val.extend([1, 1])
                variables.append(f"N_{i}_{j}_fibre")
                val.append(-(1 - f_switch) * int(self.g.edges[[i, j]]["capacity"]) * no_channels)
                lin_expressions = [cplex.SparsePair(ind=variables, val=val)]
                self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["L"], rhs=[0.])

    def source_more_than_detectors_approx(self):
        """
        Adds the constraint that approximates the requirement for more sources than detectors on every adjacent node:
        N_{j}^{D} \leq N_{i}^{S} \forall j \in N(i)
        """
        for i in self.g.nodes:
            for j in self.g.adj[i]:
                variables = [f"N_{j}_D", f"N_{i}_S"]
                val = [1, -1]
                lin_expressions = [cplex.SparsePair(ind=variables, val=val)]
                self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["L"], rhs=[0.])

    def add_fibre_on_constraint(self, Lambda):
        """
        adds the constraint that a fibre must be on to have any fibres N_{i,j}^{fibres} \leq Lambda \delta_{i,j}
        """
        for i, j in self.g.edges:
            if i < j:
                variables = [f"N_{i}_{j}_fibre", f"delta_{i}_{j}"]
                val = [1, -Lambda]
                lin_expressions = [cplex.SparsePair(ind=variables, val=val)]
                self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["L"], rhs=[0.])

    def maximise_lowest_capacity_connection_objective(self, Tij):
        """
        Adds the maximise lowest capacity objective:
        maximise z
        subject to:
        z \leq \sum_{n \in N(j)} x_{(n,j)}^{k=(i,j)} + x_{(j,n)}^{k=(j,i)}
        """
        variable_names = ['z']
        self.prob.variables.add(names=variable_names, types=[self.prob.variables.type.continuous])
        for k in self.key_dict:
            variables = []
            val = []
            for n in self.g.neighbors(k[1]):
                if Tij == None:
                    variables.extend([f"x{n}_{k[1]}_k{k[0]}_{k[1]}", f"x{k[1]}_{n}_k{k[1]}_{k[0]}"])
                    val.append(-1)
                    val.append(-1)
                elif Tij[k] > 0.000001:
                    variables.extend([f"x{n}_{k[1]}_k{k[0]}_{k[1]}", f"x{k[1]}_{n}_k{k[1]}_{k[0]}"])
                    val.append(-1 / Tij[k])
            variables.append("z")
            val.append(1)
            lin_expressions = [cplex.SparsePair(ind=variables, val=val)]
            self.prob.linear_constraints.add(lin_expr=lin_expressions, senses=["L"], rhs=[0.])
        obj_vals = [("z", 1)]
        self.prob.objective.set_linear(obj_vals)
        self.prob.objective.set_sense(self.prob.objective.sense.maximize)



    def trusted_node_optimisation_program(self, time_limit=1e5, csource=0.1, cdet=1,
                                                                c_dist=0.01, c_conn_dist=0.1, Lambda=10, cgraph=100,
                                                                Nmax=10, f_switch=0.1, no_channels = 1, Tij = None):
        """
        set up and solve the problem for minimising the overall cost of the network
        """
        print("Start Optimisation")
        self.flow_conservation_constraint()
        # add_capacity_constraint(prob, g, key_dict, Lambda=Lambda)
        self.no_flow_into_source_out_sink()
        cfibres, cconn = get_cfibre(self.g, c_dist, c_conn_dist)
        self.cost_constraint(csource, cdet, cfibres, cconn, cgraph, Nmax)
        # add_minimise_trusted_nodes_objective(prob, g)
        self.add_source_capacity_on_edge_constraint(f_switch)
        self.add_detector_capacity_on_edge_constraint(f_switch)
        self.fibre_constraint(f_switch, no_channels=no_channels)
        # self.source_more_than_detectors_approx()
        self.add_fibre_on_constraint(Lambda)
        self.maximise_lowest_capacity_connection_objective(Tij)

        self.prob.write("test_multi.lp")
        self.prob.parameters.lpmethod.set(3)
        self.prob.parameters.mip.limits.cutpasses.set(1)
        self.prob.parameters.mip.strategy.probe.set(-1)
        self.prob.parameters.mip.strategy.variableselect.set(2)
        self.prob.parameters.mip.strategy.kappastats.set(1)
        self.prob.parameters.mip.tolerances.mipgap.set(float(0.01))
        # prob.parameters.simplex.limits.iterations = 50
        print(self.prob.parameters.get_changed())
        self.prob.parameters.timelimit.set(time_limit)
        t_1 = time.time()

        self.prob.solve()
        t_2 = time.time()
        print("Time to solve problem: " + str(t_2 - t_1))
        print(f"The minimum Cost of Network: {self.prob.solution.get_objective_value()}")
        print(f"Number of Variables = {self.prob.variables.get_num()}")
        print(f"Number of Conditions = {self.prob.linear_constraints.get_num()}")
        sol_dict = self.create_sol_dict()
        flow_dict, numbers_dict, delta_dict = self.split_sol_dict(sol_dict)
        return sol_dict, self.prob



def get_cfibre(g, c_dist, c_conn_dist):
    """
    get the cost of the fibre
    """
    cfibre = {}
    cconn = {}
    for i,j in g.edges:
        if i < j:
            distance = g.edges[[i,j]]["distance"]
            cfibre[i,j] = distance * c_dist
            cconn[i,j] = distance * c_conn_dist
    return cfibre, cconn


def make_key_dict_bidirectional(key_dict):
    """
    make the key dict bidirectional : if (source, target) in key_dict then (target, source) should be too
    """
    missing_entries = [(k[1],k[0]) for k in key_dict if (k[1],k[0]) not in key_dict]
    for idx in missing_entries:
        key_dict[idx] = 0
    return key_dict#


def plot_graph(graph, delta_dict):
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
    plt.show()


if __name__ == "__main__":
    key_dict, g = utils.import_problem_trusted_nodes(node_file_path="small_node_graph_test_nodes.csv",edge_file_path = "small_node_graph_test.csv")
    for key in g.keys():
        key_dict_new = make_key_dict_bidirectional(key_dict[key])
        prob = cplex.Cplex()
        optimisation = Switched_Trusted_Node_Optimisation_Constant_Switching_Ratio(prob, g[key], key_dict_new)
        sol_dict, prob = optimisation.trusted_node_optimisation_program(Nmax =24, cgraph= 40)
        flow_dict, numbers_dict, delta_dict = optimisation.split_sol_dict(sol_dict)
        plot_graph(g[key], delta_dict)