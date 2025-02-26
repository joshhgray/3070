
def rank_based_selection(population, selection_cutoff, population_size):
    """
    Performs a rank based selection of which cutoff percent of molecules
    will be selected based on rank which is calculated based on fitness
    
    :param population:
    :param selection_cutoff:
    :returns:
    """
    mols = [node for node in population.nodes if population.nodes[node]["level"] == "Individual"]

    # Rank mols by fitness score
    ranked_mols = sorted(mols,
                         key=lambda x: population.nodes[x]["raw_fitness"],
                         reverse=True)
    
    # Select top individuals by highest fitness
    parents = ranked_mols[: len(ranked_mols) // int(100 / selection_cutoff)]

    # Prune least fit individuals (Death)
    if len(population.nodes) > population_size:
        lowest_fitness_mols = ranked_mols[len(parents):]
        for mol in lowest_fitness_mols:
            population.remove_node(mol)

    return parents