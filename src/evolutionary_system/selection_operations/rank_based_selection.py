from src.evolutionary_system.selection_operations.selection import Selection

class RankBasedSelection(Selection):

    def select(population, **kwargs):
        """
        Performs a rank based selection of which cutoff percent of molecules
        will be selected based on rank which is calculated based on fitness
        
        :param population:
        :param selection_cutoff:
        :returns:
        """
        selection_cutoff = kwargs.get("selection_cutoff")

        mols = [node for node in population.nodes if population.nodes[node]["level"] == "Individual"]

        # Rank mols by fitness score
        ranked_mols = sorted(mols,
                            key=lambda x: population.nodes[x]["raw_fitness"],
                            reverse=True)
        
        # Select top individuals by highest fitness
        parents = ranked_mols[: len(ranked_mols) * selection_cutoff // 100]

        return parents