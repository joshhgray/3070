"""
This includes a basic implementation of a Hierarchical encoding schema
for the Multi-Level Selection Genetic Algorithm (MLSGA). 

(In BGC mode)

Root Node: Represents whole population
Group Nodes: Represent clusters of BGCs
Individual Nodes: Represent individual BGCs
"""
from collections import defaultdict
import networkx as nx
from src.data_pipeline.sample_population import sample_bgcs, sample_smiles
from rdkit import Chem
import pandas as pd
import uuid
import os

def make_population_graph(population_size, population_type):
    """
    Build a hierarchical graphical representation of the population based on 
    the BGC JSON data

    :param population_size: Size of initial population.
    :returns: NetworkX DiGraph representing the hierarchical population structure.
    """
    # Iniitalize the directed graph
    population_graph = nx.DiGraph()

    # Add Root node (Population-level)
    population_graph.add_node(
        "Population",
        level="Root",
        diversity=0.0,
        total_fitness=0.0,
    )
    # SMILES population building logic
    if population_type == "smiles":
        # Read preprocessed SMILES CSV
        cd = os.path.dirname(os.path.abspath(__file__))
        preprocessed_smiles = os.path.join(cd, "../../preprocessed_smiles.csv")
        df = pd.read_csv(preprocessed_smiles)

        # Extract selected number of SMILES from the pool and convert to python list
        df = df[df["Smiles"].notna()]
        sampled_smiles = df.sample(n=population_size, replace=True)["Smiles"].tolist()

        group_name = "SMILES_group" # TODO - temp
        population_graph.add_node(group_name, level="Group", size=population_size)
        population_graph.add_edge("Population", group_name)

        for smiles in sampled_smiles:
            mol = Chem.MolFromSmiles(smiles)
            # Ensure chemical validity before adding to population
            if mol is None:
                continue
            try:
                Chem.SanitizeMol(mol)
            except:
                continue
            
            node_id = f"mol_{uuid.uuid4().hex[:10]}"
            population_graph.add_node(
                node_id,
                level="Individual",
                raw_fitness=0.0,
                taxonomy=None,
                biosynthesis_class=None,
                compounds=[{"structure": smiles, "mol": mol}],
                elite=False
            )
            population_graph.add_edge(group_name, node_id)

    # BGC population building logic
    elif population_type == "bgc":
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

        # Add Group and Inidividual-level nodes
        for group_name, group_bgcs in groups.items():
            population_graph.add_node(
                group_name,
                level="Group",
                group_fitness=0.0,
            )
            population_graph.add_edge("Population", group_name)

            for bgc in group_bgcs:
                # Extract relevant attributes from bgc data
                accession_number = bgc.get("accession")
                taxonomy = bgc.get("taxonomy", {}).get("name", "Unknown")
                compounds = bgc.get("compounds", [])

                # Format the compound(s) data
                compound_data = []
                for compound in compounds:
                    structure = compound.get("structure", "Unknown")
                    mass = compound.get("mass", 0.0)
                    formula = compound.get("formula", "Unknown")
                    
                    # Extract RDKit Mol object from SMILES string to store in graph.
                    mol = Chem.MolFromSmiles(structure) if structure != "Unknown" else None
                    # Ensure validity
                    if mol is None:
                        continue
                    try:
                        Chem.SanitizeMol(mol)
                    except:
                        continue

                    compound_data.append({
                        "structure": structure,
                        "mass": mass,
                        "formula": formula,
                        "mol": mol
                    })

                # Add the node
                population_graph.add_node(
                    accession_number,
                    level="Individual",
                    raw_fitness=0.0,
                    accession_number=accession_number,
                    biosynthesis_class=group_name,
                    taxonomy=taxonomy,
                    compounds=compound_data,
                    elite=False,
                )
                population_graph.add_edge(group_name, accession_number)

    return population_graph