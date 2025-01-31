"""
This is a basic implementation of a Hierarchical encoding schema
for the Multi-Level Selection Genetic Algorithm (MLSGA).

Root Node: Represents whole population
Group Nodes: Represent clusters of individuals
Individual Nodes: Represent individual compounds
"""
from collections import defaultdict
#from rdkit import RDLogger
import networkx as nx
import numpy as np
import json
import os

def make_population_graph(bgc_data):
    """
    Build a hierarchical graphical representation of the population based on 
    the BGC JSON data

    :param bgc_data: List of parsed BGC data.
    :returns: NetworkX DiGraph representing the hierarchical population structure
    """
    # TODO - Temporary suppression of rdlogger
    #RDLogger.DisableLog('rdApp.*')

    # Group BGCs by biosynthesis class
    # Default missing attributes to "Unknown"
    groups = defaultdict(list)
    for bgc in bgc_data:
        biosynthesis_classes = [
            biosynthesis_class.get("class", "Unknown") 
            for biosynthesis_class in bgc.get("biosynthesis", {}).get("classes", [])
        ]
        if not biosynthesis_classes:
            biosynthesis_classes = ["Unkown"]
        for biosynthesis_class in biosynthesis_classes:
            groups[biosynthesis_class].append(bgc)

    # Iniitalize the directed graph
    population_graph = nx.DiGraph()

    # Add Root node (Population-level)
    population_graph.add_node(
        "Population",
        level="Root",
        diversity=0.0,
        total_fitness=0.0,
        total_mass=0.0,
        total_cyclic_count=0 # TODO - decide if relevant
    )

    # Add Group and Inidividual-level nodes
    """
    TODO - I've added some extra attributes here, I don't have a plan to use
    them but if there is remote possibility I could use them I will keep them
    so I can use them later on.
    """
    for group_name, group_bgcs in groups.items():
        population_graph.add_node(
            group_name,
            level="Group",
            total_fitness=0.0,
            average_fitness= 0.0,
            size=len(group_bgcs),
            total_mass=0.0,
            average_mass=0.0,
            total_cyclic_count=0,
        )
        population_graph.add_edge("Population", group_name)

        for bgc in group_bgcs:
            # Extract relevant attributes from bgc data
            accession_number = bgc.get("accession")
            biosynthesis_class = [
                biosynthesis_class.get("class", "Unknown")
                for biosynthesis_class in bgc.get("biosynthesis", {}).get("classes", [])
            ]
            structure = bgc.get("compounds", [{}])[0].get("structure", "Unknown")
            mass = bgc.get("compounds", [{}])[0].get("mass", 0.0)
            formula = bgc.get("compounds", [{}])[0].get("formula", "Unknown")

            taxonomy = bgc.get("taxonomy", {}).get("name", "Unknown")

            # Add the node
            population_graph.add_node(
                accession_number,
                level="Individual",
                qed=0.0,
                raw_fitness=0.0,
                scaled_fitness=0.0,
                accession_number=accession_number,
                biosynthesis_class=biosynthesis_class,
                structure=structure,
                mass=mass,
                formula=formula,
                taxonomy=taxonomy
            )
            population_graph.add_edge(group_name, accession_number)
    return population_graph










