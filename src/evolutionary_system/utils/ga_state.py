# Centralized State Management of Genetic Algorithm
ga_active = False
latest_population = None
latest_diversity_log = []

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