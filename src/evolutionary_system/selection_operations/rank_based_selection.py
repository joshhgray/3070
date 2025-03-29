from src.evolutionary_system.selection_operations.selection import Selection

class RankBasedSelection(Selection):

    def select(self, population, **kwargs):
        """
        Performs a rank based selection on a given set of nodes (a group).
        Defaults to full population if no group is provided.

        :param population: NetworkX graph representign the population.
        :param selection_cutoff: percentage cutoff to select.
        :parm nodes (optional): Group set of individuals (their IDs).
        :returns: List of selected parents node IDs
        """
        selection_cutoff = kwargs.get("selection_cutoff")
        nodes = kwargs.get("nodes")

        # Full Population
        if nodes is None:
            mols = [node for node in population.nodes if population.nodes[node].get("level") == "Individual"]
        # Group
        else:
            mols = [node for node in nodes if population.nodes[node].get("level") == "Individual"]

        # Rank mols by fitness score
        ranked_mols = sorted(mols,
                            key=lambda node: population.nodes[node].get("raw_fitness"),
                            reverse=True)
        
        # Select only top individuals 
        cutoff = max(1, len(ranked_mols) * selection_cutoff // 100)
        return ranked_mols[:cutoff]