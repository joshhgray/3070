"""
This is a basic implementation of a Hierarchical encoding schema
for the Multi-Level Selection Genetic Algorithm (MLSGA).

Root Node: Represents whole population
Group Nodes: Represent clusters of individuals
Individual Nodes: Represent individual compounds
"""

import networkx as nx # type: ignore

"""
Input: population - list of SMILES strings
Output: Graph representing hierarchy of population - an nx.DiGraph
"""
def make_graph(population):

    # Total number of groups
    num_groups = 4
    individuals_per_group = len(population) // num_groups

    # Initialize graph
    population_graph = nx.DiGraph()

    # Initialize root node
    population_graph.add_node("Population", level="Root", diversity=None)

    # Initialize group and individual nodes
    for group_id in range(1, num_groups+1):
        # Initalize group nodes
        group_name = f"Group_{group_id}"
        population_graph.add_node(group_name, 
                        level="Group", 
                        total_fitness=None, 
                        average_fitness=None) 
        population_graph.add_edge("Population", group_name)

        # Create evenly sized groups
        group_start = (group_id - 1) * individuals_per_group
        group_end = group_start + individuals_per_group
        group = population[group_start:group_end]

        # Initizlize individual nodes
        for individual_id, smiles_str in enumerate(group, start=1):
            individual_name = f"Individual_{group_id}_{individual_id}"
            population_graph.add_node(individual_name,
                            level="Individual",
                            smiles_str=smiles_str,
                            qed=None,
                            raw_fitness=None,
                            scaled_fitness=None)
            population_graph.add_edge(group_name, individual_name)

    return population_graph