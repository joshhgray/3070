# Centralized State Management of Genetic Algorithm
ga_active = False
latest_population = None
latest_diversity_log = []
mutation_success_rate_log = []
crossover_success_rate_log = []
latest_mutation_success_rates = 0
latest_crossover_success_rates = 0

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