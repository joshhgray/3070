from abc import ABC, abstractmethod

class Crossover(ABC):
    """
    This is the base crossover class for all types of crossover operations.
    """

    @abstractmethod
    def apply(self, parent1, parent2):
        """Apply crossover to a set of two RDKit Mol Objects and return offspring."""
        pass