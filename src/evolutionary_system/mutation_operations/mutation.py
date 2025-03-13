from abc import ABC, abstractmethod
from rdkit import Chem

class Mutation(ABC):
    """
    This is the base mutaiton class for all types of mutation operations.
    """

    @abstractmethod
    def apply(self, rw_mol):
        """Apply mutation to an RDKit RWMol and return RDKit Mol"""
        pass