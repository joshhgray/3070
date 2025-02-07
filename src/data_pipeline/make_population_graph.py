"""
This is a basic implementation of a Hierarchical encoding schema
for the Multi-Level Selection Genetic Algorithm (MLSGA).

Root Node: Represents whole population
Group Nodes: Represent clusters of BGCs
Individual Nodes: Represent individual BGCs
"""
from collections import defaultdict
#from rdkit import RDLogger
import networkx as nx
from src.data_pipeline.sample_population import sample_bgcs

def make_population_graph(population_size):
    """
    Build a hierarchical graphical representation of the population based on 
    the BGC JSON data

    :param population_size: Size of initial population.
    :returns: NetworkX DiGraph representing the hierarchical population structure.
    """

    sampled_bgcs = sample_bgcs(population_size)

    # Group BGCs by biosynthesis class
    # Default missing attributes to "Unknown"
    groups = defaultdict(list)
    for bgc in sampled_bgcs:
        biosynthesis_classes = [
            biosynthesis_class.get("class", "Unknown") 
            for biosynthesis_class in bgc.get("biosynthesis", {}).get("classes", [])
        ]

        if not biosynthesis_classes:
            biosynthesis_classes = ["Unknown"]
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
    )

    # Add Group and Inidividual-level nodes
    """
    TODO - I've added some extra attributes here, I don't have a plan to use
    them but if there is remote possibility I could use them I will keep them
    so I can use them later on. Attributes can be commented out.
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
            taxonomy = bgc.get("taxonomy", {}).get("name", "Unknown")
            compounds = bgc.get("compounds", [])

            # Format the compound(s) data
            compound_data = []
            total_mass = 0.0
            for compound in compounds:
                structure = compound.get("structure", "Unkonwn")
                mass = compound.get("mass", 0.0)
                formula = compound.get("formula", "Unkonw")
                mol_graph = compound.get("mol_graph", None)

                compound_data.append({
                    "structure": structure,
                    "mass": mass,
                    "formula": formula,
                    "mol_graph": mol_graph
                })

                total_mass += mass
            average_mass = total_mass / len(compound_data) if compound_data else 0.0

        
            # Add the node
            population_graph.add_node(
                accession_number,
                level="Individual",
                qed=0.0,
                raw_fitness=0.0,
                scaled_fitness=0.0,
                accession_number=accession_number,
                biosynthesis_class=biosynthesis_class,
                taxonomy=taxonomy,
                compounds=compound_data,
                average_mass=average_mass
            )
            population_graph.add_edge(group_name, accession_number)
    # si hay un error cuando haciendo el graph en el initial population
    # for node in population_graph.nodes:
    #     if population_graph.nodes[node]["level"] == "Individual":
    #         for compound in population_graph.nodes[node]["compounds"]:
    #             if "mol_graph" not in compound or compound["mol_graph"] is None:
    #                 print(f"Missing mol_graph for {node}")
    return population_graph










