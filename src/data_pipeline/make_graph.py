"""
This is a basic implementation of a Hierarchical encoding schema
for the Multi-Level Selection Genetic Algorithm (MLSGA).

Root Node: Represents whole population
Group Nodes: Represent clusters of individuals
Individual Nodes: Represent individual compounds
"""
from rdkit import RDLogger # type: ignore
import networkx as nx
import numpy as np # type: ignore
"""
Input: population - list of SMILES strings
Output: Graph representing hierarchy of population - an nx.DiGraph
"""
def make_graph(population):
    # TODO - Temporary suppression of rdlogger
    RDLogger.DisableLog('rdApp.*')

    # Total number of groups
    num_groups = 5 # TODO - make dynamic
    groups = np.array_split(population, num_groups)

    # Initialize graph
    population_graph = nx.DiGraph()

    # Initialize root node
    population_graph.add_node("Population", level="Root", diversity=0.0)

    # Initialize group and individual nodes
    for group_id, group in enumerate(groups, start=1):
        # Convert back to list - not sure if neccessary
        group = list(group)
        
        # Initalize group nodes
        group_name = f"Group_{group_id}"
        population_graph.add_node(group_name, 
                        level="Group", 
                        total_fitness=0.0, 
                        average_fitness=0.0) 
        population_graph.add_edge("Population", group_name)

        # Initizlize individual nodes
        for individual_id, smiles_str in enumerate(group, start=1):
            individual_name = f"Individual_{group_id}_{individual_id}"
            population_graph.add_node(individual_name,
                            level="Individual",
                            smiles_str=smiles_str,
                            qed=0.0,
                            raw_fitness=0.0,
                            scaled_fitness=0.0)
            population_graph.add_edge(group_name, individual_name)

    return population_graph