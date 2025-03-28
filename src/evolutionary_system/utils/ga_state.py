# Centralized State Management of Genetic Algorithm
ga_active = False
latest_population = None
latest_diversity_log = []
mutation_success_rate_log = []
crossover_success_rate_log = []
latest_mutation_success_rates = 0
latest_crossover_success_rates = 0
current_population_size = 0
current_generation_number = 0
mol_weight_threshold = 500 # TODO - move to config
selected_fitness_functions = [] # TODO - can be infered from config rather than in here
fitness_history_log = []
active_filters = []
tuning_weights = {}
mw_range = [250, 450]
logp_range = [0.5, 5.0]

MUTATION_PROBABILITIES = {
    "atomic_substitution": 0.25,
    "functional_group": 0.25,
    "ring": 0.25,
    "fragment": 0.25,
}

def set_ga_active(value):
    global ga_active
    ga_active = value

def is_ga_active():
    return ga_active

def update_latest_population(population):
    global latest_population
    latest_population = population

def get_latest_population():
    return latest_population

def update_latest_diversity(diversity_log):
    global latest_diversity_log
    latest_diversity_log = diversity_log

def get_latest_diversity():
    return latest_diversity_log

def update_latest_mutation_rates(rates):
    global latest_mutation_success_rates
    latest_mutation_success_rates = rates

def get_latest_mutation_rates():
    return latest_mutation_success_rates

def update_latest_crossover_rates(rates):
    global latest_crossover_success_rates
    latest_crossover_success_rates = rates

def get_latest_crossover_rates():
    return latest_crossover_success_rates

def update_mutation_log(rate):
    global mutation_success_rate_log
    mutation_success_rate_log.append(rate)

def get_mutation_log():
    return mutation_success_rate_log

def update_crossover_log(rate):
    global crossover_success_rate_log
    crossover_success_rate_log.append(rate)

def get_crossover_log():
    return crossover_success_rate_log

def update_current_population_size(population_size):
    global current_population_size
    current_population_size = population_size

def get_current_population_size():
    return current_population_size

def update_current_generation_number(generation_number):
    global current_generation_number
    current_generation_number = generation_number

def get_current_generation_number():
    return current_generation_number

def update_mol_weight_threshold(value):
    global mol_weight_threshold
    mol_weight_threshold = value

def get_mol_weight_threshold():
    return mol_weight_threshold

def update_selected_fitness_functions(functions):
    # Takes list of fitness functions
    global selected_fitness_functions
    selected_fitness_functions = functions

def get_selected_fitness_functions():
    return selected_fitness_functions

def update_mutation_probabilities(updated_probabilities):
    global MUTATION_PROBABILITIES
    MUTATION_PROBABILITIES = updated_probabilities

def get_mutation_probabilities():
    return MUTATION_PROBABILITIES

def update_fitness_log(fitness_tuple):
    global fitness_history_log
    fitness_history_log.append(fitness_tuple)

def get_fitness_history_log():
    return fitness_history_log

def update_active_filters(filters):
    global active_filters
    active_filters = filters

def get_active_filters():
    return active_filters

def update_tuning_weights(weights):
    global tuning_weights
    tuning_weights = weights

def get_tuning_weights():
    return tuning_weights

def update_mw_range(range):
    global mw_range
    mw_range=range

def get_mw_range():
    return mw_range

def get_mw_target():
    # Define target as mid point on range
    low, high = get_mw_range()
    return (low + high) / 2

def update_logp_range(range):
    global logp_range
    logp_range=range

def get_logp_range():
    return logp_range
