from abc import ABC, abstractmethod

class Selection(ABC):
    """
    This is the base selection class for all types of selection methods.
    """

    @abstractmethod
    def select(self, population, **kwargs):
        """Selects a subset of individual node IDs from the population for breeding."""
        pass