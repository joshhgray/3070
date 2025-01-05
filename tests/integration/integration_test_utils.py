import pytest #type:ignore

@pytest.fixture(scope="session")
def mock_config():
    return {
        "crossover_rate": 0.1,
        "fitness_weights": [1.0, 1.0],
        "mutation_rate": 0.1,
        "num_elite_groups": 1,
        "num_elite_individuals": 1,
        "num_generations": 10,
        "num_threads": 1,
        "population_size": 50,
        "selection_method": "sus",
    }

@pytest.fixture(scope="session")
def mock_population():
    return [
        "c1(Cl)c(c[nH]c1Cl)[N+](=O)[O-]",
        "c1c(c(c([nH]1)Cl)Cl)[N+](=O)[O-]",
        "Clc1c([N+]([O-])=O)c[nH]c1Cl",
        "[O-][N+](c1c(c(Cl)[nH]c1)Cl)=O",
        "c1(c([N+](=O)[O-])c[nH]c1Cl)Cl",
        "[N+](c1c[nH]c(Cl)c1Cl)([O-])=O",
        "c1([N+](=O)[O-])c(Cl)c(Cl)[nH]c1",
        "c1(Cl)c([nH]cc1[N+]([O-])=O)Cl",
        "c1[nH]c(Cl)c(Cl)c1[N+]([O-])=O",
        "O=[N+]([O-])c1c[nH]c(Cl)c1Cl",
        "c1(c(Cl)[nH]cc1[N+]([O-])=O)Cl",
        "c1(Cl)[nH]cc([N+]([O-])=O)c1Cl",
        "c1(Cl)c(Cl)c(c[nH]1)[N+](=O)[O-]",
        "c1([N+](=O)[O-])c(c(Cl)[nH]c1)Cl",
        "Clc1[nH]cc(c1Cl)[N+]([O-])=O",
        "Clc1c(Cl)[nH]cc1[N+](=O)[O-]",
        "c1c(c(c(Cl)[nH]1)Cl)[N+]([O-])=O",
        "Clc1c(c([N+]([O-])=O)c[nH]1)Cl",
        "c1(c(c(Cl)[nH]c1)Cl)[N+]([O-])=O",
        "c1c(c(Cl)c(Cl)[nH]1)[N+]([O-])=O",
    ]