import pytest #type:ignore
from src.evolutionary_system.ga import run_ga
from src.data_pipeline.make_graph import make_graph

def test_run_ga(mock_population, mock_config):
    '''
    Integration testing for the Genetic Algorithm (GA)
    '''
    # Setup
    test_population = make_graph(mock_population)

    final_population, diversity_log = run_ga(
        initial_population=test_population,
        population_size=mock_config["population_size"],
        num_generations=mock_config["num_generations"],
        mutation_rate=mock_config["mutation_rate"],
        crossover_rate=mock_config["crossover_rate"],
        num_elite_individuals=mock_config["num_elite_individuals"],
        num_elite_groups=mock_config["num_elite_groups"],
        selection_method=mock_config["selection_method"],
        fitness_weights=mock_config["fitness_weights"],
        num_threads=mock_config["num_threads"]
    )

    # Tests

    # TODO - test final pop is graph
    # TODO - test fitness eval
    # TODO - test diversity log
    # TODO - test parameters

    # TODO - temp
    assert final_population is not None