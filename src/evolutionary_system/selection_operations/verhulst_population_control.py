def verhulst_population_control(population_size, carrying_capacity):
    """
    Control population size using Verhulst's model of population growth
    dN/dt = rN(1-N/K)
    
    :param population_size: Total number of molecules in population (N)
    :param carrying_capacity: Max population size (K) - 
    TODO - need to find a K that isn't too computationally expensive but still allows exploratory growth 
    :returns: float from 0-1 representing survival probability for mols in current generation
              (aka scale of population growth for current generation)
    """
    # TODO - find good scales - trying out an explicit limit
    if population_size >= carrying_capacity:
        return 0
    else:
        # Minimum 0.4 (40%) survival probability
        return max(0.4, 1 - (population_size / (1.5 * carrying_capacity)))