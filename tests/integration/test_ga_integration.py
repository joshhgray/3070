import pytest #type:ignore
from src.evolutionary_system.ga import run_ga
from src.data_pipeline.make_population_graph import make_graph

# TODO - break up into separate functions
def test_ga(mock_population, mock_config):
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

    # Validation of outputs
    # Final Population
    assert final_population is not None

    individuals = [
        node for node in final_population.nodes
        if final_population.nodes[node]["level"] == "Individual"
    ]

    # Test QED calculation
    for individual in individuals:
        qed_score = final_population.nodes[individual].get("qed")
        assert qed_score is not None
        assert isinstance(qed_score, float)
        assert 0.0 <= qed_score <= 1.0

    # Diversity Log
    # Test Individual Level Diversity Calculation
    assert diversity_log is not None
    assert len(diversity_log) == mock_config["num_generations"]
    assert all(isinstance(diversity_score, float) for diversity_score in diversity_log)
    assert all(0.0 <= diversity_score <= 1.0 for diversity_score in diversity_log)
    
    # Test Population Level Diversity Calculation
    population_diversity = final_population.nodes["Population"].get("diversity")
    assert population_diversity is not None
    assert isinstance(population_diversity, float)
    assert 0.0 <= population_diversity <= 1.0