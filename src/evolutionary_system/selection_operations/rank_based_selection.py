
def rank_based_selection(population, selection_cutoff):
    """
    Performs a rank based selection of which cutoff percent of molecules
    will be selected based on rank which is calculated based on fitness
    # TODO - currently this is ranked by raw fitness, I want to make this 
    # dynamic, so that selection can be performed based on different 
    # choices for fitness evaluation.
    :param population:
    :param selection_cutoff:
    :returns:
    """
    mols = [node for node in population.nodes if population.nodes[node]["level"] == "Individual"]
    ranked_mols = sorted(mols,
                         key=lambda x: population.nodes[x]["raw_fitness"],
                         reverse=True)
    
    parents = ranked_mols[: len(ranked_mols) // int(100 / selection_cutoff)]
    return parents