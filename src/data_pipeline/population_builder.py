import os
import json
import random
import networkx as nx

def load_bgc_jsons(json_dir):
    """
    
    """
    bgc_data = []
    for file in os.listdir(json_dir):
        if file.endswith(".json"):
            file_path = os.path.jion(json_dir, file)
            with open(file_path, "r") as f:
                data = json.load(f)
                bgc_data.append(data)

    return bgc_data

def build_bgc_graph(bgc):
    """
    
    """
    graph = nx.Graph()

    # Cluster node
    cluster_id = bgc.get("cluster_id")
    graph.add_node(cluster_id, type="cluster", description=bgc.get("description"))

    # Gene Nodes
    for i, gene in enumerate(bgc.get("genes")):
        gene_id = f"{cluster_id}_gene_{i}"
        gene_sequence = gene.get("sequence")
        gene_function = gene.get("product")

        # Connect to cluster
        graph.add_node(gene_id, type="gene", sequence=gene_sequence, function=gene_function)
        graph.add_edge(cluster_id, gene_id)
    return graph

def build_population(json_dir, population_size):
    """
    """
    bgc_data = load_bgc_jsons(json_dir)
    # randomly sample population from pool
    # TODO - make this seedable
    sampled_bgcs = random.sample(bgc_data, min(population_size, len(bgc_data)))

    population = [build_bgc_graph(bgc) for bgc in sampled_bgcs]
    return population

if __name__ == "__main__":
    json_dir = "../mibig_json_4.0"