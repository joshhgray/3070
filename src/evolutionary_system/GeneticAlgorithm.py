from src.evolutionary_system.fitness_operations.aggregate_fitness import aggregate_fitness
from src.evolutionary_system.fitness_operations.calculate_population_diversity import calculate_population_diversity
from src.evolutionary_system.selection_operations.rank_based_selection import RankBasedSelection
from src.evolutionary_system.selection_operations.stochastic_universal_sampling import StochasticUniversalSampling
from src.evolutionary_system.selection_operations.verhulst_population_control import verhulst_population_control
from src.evolutionary_system.mutation_operations.bioisosteric_mutations import BioisostericMutation
from src.evolutionary_system.mutation_operations.ring_mutation import RingMutation
from src.evolutionary_system.mutation_operations.fragment_mutaiton import FragmentMutation
from src.evolutionary_system.crossover_operations.max_common_substruct_crossover import MCSCrossover
from src.evolutionary_system.mutation_operations.atomic_mutations.atomic_substitution import AtomicSubstitutionMutation
from src.data_pipeline.make_population_graph import make_population_graph
from src.evolutionary_system.utils.config_loader import load_config
from src.evolutionary_system.utils.logger import log_metrics
from src.evolutionary_system.utils.ga_state import (
    update_latest_diversity, update_latest_population, is_ga_active, 
    set_ga_active, update_latest_crossover_rates, update_latest_mutation_rates, 
    update_crossover_log, update_mutation_log, update_current_population_size, 
    update_current_generation_number, update_fitness_log, get_tuning_weights,
    get_active_filters, get_mw_range, get_logp_range, get_fitness_history_log,
    get_mutation_log, get_crossover_log
)
from rdkit import Chem
import pandas as pd
import numpy as np
import datetime
import psutil
import random
import uuid
import time
import csv
import os

CROSSOVER_METHODS = {
    "mcs_crossover": MCSCrossover(),
}
SELECTION_METHODS = {
    "rank_based_selection": RankBasedSelection(),
    "stochastic_universal_sampling": StochasticUniversalSampling(),
}
FITNESS_CONFIG = {
    "tuning_weights": get_tuning_weights(),
    "active_filters": get_active_filters(),
    "mw_range": get_mw_range(),
    "logp_range": get_logp_range()
}

class GeneticAlgorithm:
    def __init__(self, config):
        self.config = config
        self.population_size = config["population_size"]
        self.num_generations = config["num_generations"]
        self.crossover_methods = config.get("crossover_methods")
        self.elitism_weight = config.get("elitism_weight")
        self.mutation_weights = config.get("mutation_weights")
        self.selection_method = config.get("selection_method")
        self.selection_cutoff = 0.5
        self.carrying_capacity = config["carrying_capacity"]
        self.population_type = config["data_source"]
        self.working_population = None
        self.diversity_log = []
        self.fragment_library = self.load_fragment_library()
        self.fragment_mutator = FragmentMutation(self.fragment_library)

    def load_fragment_library(self):
        """
        # Loads the fragment library (CSV) and saves fragment's SMILES to a list.
        """
        cd = os.path.dirname(os.path.abspath(__file__))
        lirbrary_path = os.path.join(cd, "../../bgc_fragment_library.csv")
        df = pd.read_csv(lirbrary_path)

        frag_list = []
        for _, row in df.iterrows():
            frag_smiles = row["fragment_smiles"]
            frag_list.append(frag_smiles)
        
        return frag_list

    def initialize_population(self):
        """
        # TODO
        """
        self.working_population = make_population_graph(self.population_size, self.population_type)
        self.evaluate_initial_fitness()

    def evaluate_initial_fitness(self):
        """
        # TODO
        """
        for node in self.working_population.nodes:
            if self.working_population.nodes[node]["level"] == "Individual":
                compounds = self.working_population.nodes[node].get("compounds")
                mol = compounds[0].get("mol") if compounds else None
                self.working_population.nodes[node]["raw_fitness"] = aggregate_fitness(mol)
    
    def apply_selection(self):
        """
        Multi-level selection where groups are selected via SUS based on group fitness,
        then apply standard selection within selected groups
        """
        """
        Compute normalized, average group fitness for each group in the population.
        """
        groups = list(self.working_population.successors("Population"))
        group_fitnesses = []
        for group in groups:
            individuals = list(self.working_population.successors(group))
            # Can't seem to get NoneType's from getting into individual_fitnesses
            # So filtering out everything except ints and floats here
            individual_fitnesses = []
            for node in individuals:
                fitness_score = self.working_population.nodes[node].get("raw_fitness")
                if isinstance(fitness_score, (int, float)): # Ensure int or float type
                    individual_fitnesses.append(fitness_score)

            if individual_fitnesses:
                avg_fitness = np.mean(individual_fitnesses)
            else:
                avg_fitness = 0.0
            # Add to graph
            self.working_population.nodes[group]["group_fitness"] = avg_fitness
            group_fitnesses.append(avg_fitness)

        # Normalize group fitness
        total_fitness = sum(group_fitnesses)
        if total_fitness == 0: # avoid divide zero
            probabilities = [1 / len(groups)] * len(groups)
        else:
            probabilities = [fit_score / total_fitness for fit_score in group_fitnesses]

        """
        Group-level Selection - select groups

        Based off of: https://en.wikipedia.org/wiki/Stochastic_universal_sampling
        Original Source: 
        Baker, J. E. (1987, July). Reducing bias and inefficiency in the selection algorithm. 
        In Proceedings of the second international conference on genetic algorithms (Vol. 206,
        pp. 14-21).
        """
        selected_groups = []
        # minimum of 1 to avoid zero division
        num_selected_groups = max(1, int(len(groups) * self.selection_cutoff))

        # Set up SUS pointers
        p = 1.0 / num_selected_groups 
        start = random.uniform(0, p)
        pointers = [start + i * p for i in range(num_selected_groups)]
        cumulative_probabilities = np.cumsum(probabilities)

        i = 0
        # Roulette Wheel Selection
        for pointer in pointers:
            while i < len(cumulative_probabilities) and cumulative_probabilities[i] < pointer:
                i += 1
            if i < len(groups): # avoid out of bounds
                selected_groups.append(groups[i])

        """
        Individual-level Selection - select individuals within selected groups
        """
        selected_parents = []
        selection_operator = SELECTION_METHODS[self.selection_method
                                               ]
        for group in selected_groups:
            individuals = list(self.working_population.successors(group))
            # Apply Rank-based
            if self.selection_method == "rank_based_selection":
                parents = selection_operator.select(self.working_population,
                                                     selection_cutoff=self.selection_cutoff,
                                                     nodes=individuals)
            # Apply SUS 
            elif self.selection_method == "stochastic_universal_sampling":
                parents = selection_operator.select(self.working_population, 
                                                    root_node=group,
                                                    selection_cutoff=self.selection_cutoff)
            else:
                return None
            
            selected_parents.extend(parents)
        return selected_parents

    def select_mutation_type(self):
        """
        Given user probabilities of mutation functions, normalize, and then select one probabilistically.
        """
        mutation_probabilities = self.mutation_weights

        # No mutation operations selected
        if not mutation_probabilities:
            return None
        
        # Normalize weights
        sum_probs = sum(mutation_probabilities.values())
        if sum_probs > 0:
            normalized_weights = {key: val / sum_probs for key, val in mutation_probabilities.items()}
        else:
            return None

        # Extract mutation types and probabilites
        types = list(normalized_weights.keys())
        probs = list(normalized_weights.values())

        return random.choices(types, weights=probs, k=1)[0]

    def apply_mutation_by_type(self, mol, mutation_type):
        """
        Applies user selected mutation type
        """
        # ATOMIC
        if mutation_type == "atomic_substitution":
            return AtomicSubstitutionMutation().apply(mol)
        # FUNCTIONAL
        elif mutation_type == "functional_group":
            return BioisostericMutation().apply(mol)

        # RING
        elif mutation_type == "ring":
            return RingMutation().apply(mol)

        # FRAGMENT
        elif mutation_type == "fragment":
            return self.fragment_mutator.apply(mol)

    def apply_mutation(self, parents):
        """
        Applies one of the selected mutations to the parent molecules.
        Continusously tracks success rate of mutation operations.
        If mutation fails, fallback to the next most probable mutation operations
        based on their weights, otherwise return the original.
        """
        mutation_attempts = 0
        mutation_successes = 0

        for parent in parents:
            # Extract mol graph from parent
            # Note: t is mportant to add backup logic becuase it is easy to 
            # create invalid compounds trhought the getnic operations
            compounds = self.working_population.nodes[parent].get("compounds")
            mol = compounds[0].get("mol") if compounds else None

            # Move on if mol is invalid
            if mol is None:
                continue
            
            selected_mutation = self.select_mutation_type()
            if selected_mutation is None:
                continue
            
            mutation_attempts += 1
            # Sort mutation types by probability of occuring (based on user given weights auto or equal weight = random selection)
            sorted_mutation_types = sorted(self.mutation_weights.items(), key=lambda item: item[1], reverse=True)
            fallback_mutation_types = [mut for mut, _ in sorted_mutation_types if mut != selected_mutation]
            attempts = [selected_mutation] + fallback_mutation_types

            for mutation_type in attempts:

                try:
                    # Get an editable version of the mol and apply selected mutation type
                    rw_mol = Chem.RWMol(mol)
                    mutated_mol = self.apply_mutation_by_type(rw_mol, mutation_type)

                    if mutated_mol and mutated_mol.GetNumAtoms() > 0:
                        # Check if molecule is framented and only keep largest if so.
                        frags = Chem.GetMolFrags(mutated_mol, asMols=True, sanitizeFrags=True)

                        if len(frags) > 1:
                            mutated_mol = max(frags, key=lambda mol: mol.GetNumAtoms())

                        Chem.SanitizeMol(mutated_mol)
                        mutated_smiles = Chem.MolToSmiles(mutated_mol)

                        # Ensure the molecule has actually mutated and is not just the same
                        parent_smiles = Chem.MolToSmiles(mol)
                        if mutated_smiles == parent_smiles:
                            continue

                        if mutated_mol:
                            new_id = f"offspring_{uuid.uuid4().hex[:10]}"
                            self.working_population.add_node(new_id, level="Individual", compounds=[])
                            self.working_population.nodes[new_id]["compounds"] = [{"structure": mutated_smiles, "mol": mutated_mol}]
                            self.working_population.nodes[new_id]["raw_fitness"] = aggregate_fitness(mutated_mol)
                            mutation_successes += 1
                            break
                except: # Mutation failed, try again if any remain.
                    continue
            
        return mutation_attempts, mutation_successes

    def apply_crossover(self, parents):
        """
        Pairs parents together and iteratively attempts to perform genetic crossover, creating new
        offspring. If validation of molecule is successful, the offspring is added to the population.
        Continuosly tracks success rate of crossover mutations.
        """
        crossover_attempts = 0
        crossover_successes = 0

        # TODO - modify once grouping is properly implemeted
        random.shuffle(parents)
        # Create pairs of parents
        parent_pairs = [(parents[i], parents[i+1]) for i in range(0, len(parents) - 1, 2)]


        for parent1, parent2 in parent_pairs:
            # Extract each parent
            # TODO - empty list fallback ensures graceful failure, could avoid this by 
            #        doing a check after every generation to purge mols with missing data.
            #        * see below todo. Patching this together to get it running for now.
            mol_1_compounds = self.working_population.nodes[parent1].get("compounds", [])
            mol_2_compounds = self.working_population.nodes[parent2].get("compounds", [])

            # TODO - * see above todo, remove once implemented. 
            if not mol_1_compounds or not mol_2_compounds:
                continue

            mol_1 = mol_1_compounds[0].get("mol") if mol_1_compounds else None
            mol_2 = mol_2_compounds[0].get("mol") if mol_2_compounds else None

            # Getting stuck on large mols while running multiprocessing
            # temp fix to filter out performing crossover on too big of mols
            if mol_1.GetNumAtoms() > 100 or mol_2.GetNumAtoms() > 100:
                continue

            crossover_attempts += 1
            
            # Currently there is only one crossover method, but this will later be a probabilistic
            # choice once others are implemented - similar to mutation is now.
            if self.crossover_methods:
                crossover_method = random.choice(self.crossover_methods)
            else:
                continue

            crossover_function = CROSSOVER_METHODS.get(crossover_method)

            if crossover_function:
                offspring = crossover_function.apply(mol_1, mol_2)
                if not offspring:
                    continue
            else:
                return crossover_attempts, crossover_successes
            
            # Confirm chemical validity, skip if offspring is invalid
            try:
                Chem.SanitizeMol(offspring)
                offspring_smiles = Chem.MolToSmiles(offspring)

                # Store new molecule in population
                new_id = f"offspring_{uuid.uuid4().hex[:10]}"
                self.working_population.add_node(new_id, level="Individual", compounds=[])
                self.working_population.nodes[new_id]["compounds"] = [{"structure": offspring_smiles, "mol": offspring}]
                self.working_population.nodes[new_id]["raw_fitness"] = aggregate_fitness(offspring)

                crossover_successes += 1
            except Exception as e:
                #print(f"Error adding new mol to population graph after applying crossover: {e}")
                continue

        return crossover_attempts, crossover_successes

    def prune_population(self):
        """
        Population control logic - probabilistic pruning of low-fitness individuals
        """
        # Get all non-elite individuals
        individuals = [node for node in self.working_population.nodes 
                if self.working_population.nodes[node]["level"] == "Individual"
                and not self.working_population.nodes[node].get("elite", False)]

        # Calculate survival probabilities using Verhlust-inspired model - capped at carrying capacity.
        survival_probability = verhulst_population_control(len(individuals), self.carrying_capacity)
        fitness_values = np.array([self.working_population.nodes[node]["raw_fitness"] for node in individuals])

        fitness_sum = fitness_values.sum()
        if fitness_sum == 0 or np.isnan(fitness_sum): # avoiding div errors
            # if fitness is unavailable weights are all equalized to 1
            normalized_fit_values = np.ones_like(fitness_values) / len(fitness_values)
        else:
            normalized_fit_values = fitness_values / fitness_sum

        # Calculate probabilities of removal
        removal_probabilities = (1 - normalized_fit_values) * (1 - survival_probability)

        individuals_to_remove = []
        for idx, node in enumerate(individuals):
            if np.random.rand() < removal_probabilities[idx]:
                individuals_to_remove.append(node)

        # Prevent over-pruning and convergence 
        num_keepers = len(individuals) - len(individuals_to_remove)
        min_population_size = int(0.5 * self.population_size)
        if num_keepers < min_population_size:
            # sort those staged for purge by probability of survival (low-to-high)
            sorted_probs = sorted(zip(removal_probabilities, individuals), key=lambda prob: prob[0], reverse=True)
            # Remove at most until population size reaches minimum bound.
            individuals_to_remove = [node for _, node in sorted_probs[:len(individuals) - min_population_size]]

        # Remove selected nodes from population graph
        self.working_population.remove_nodes_from(individuals_to_remove)

    def run_ga(self):
        # Execute the main generational loop of the GA
        #print("GA: Evaluating Initial Fitness.")
        self.evaluate_initial_fitness()

        # Generational Loop
        for generation in range(self.num_generations):
            # Confirm operational status
            if not is_ga_active():
                break

            update_current_generation_number(generation + 1)

            """
            Calculate Group Fitness
            """
            for group in self.working_population.successors("Population"):
                individuals = list(self.working_population.successors(group))
                fitnesses = [self.working_population.nodes[node]["raw_fitness"] for node in individuals]
                group_fitness = sum(fitnesses) / len(fitnesses) if fitnesses else 0.0 # prevent divide zero
                self.working_population.nodes[group]["group_fitness"] = group_fitness

            """
            Elitism
            """
            individuals = [node for node in self.working_population.nodes if self.working_population.nodes[node]["level"] == "Individual"]
            ranked = sorted(individuals, key=lambda node: self.working_population.nodes[node]["raw_fitness"], reverse=True)
            num_elites = int(len(ranked) * (self.elitism_weight / 100))

            # Remove previous generation's elite status
            for node in self.working_population.nodes:
                if self.working_population.nodes[node]["level"] == "Individual":
                    self.working_population.nodes[node]["elite"] = False
            
            # Assign new elites
            elites = ranked[:num_elites]
            for node in elites:
                self.working_population.nodes[node]["elite"] = True

            """
            Genetic Operations
            """
            # Selection
            #print("DEBUG: Selecting Parents")
            parents = self.apply_selection()

            #print("DEBUG: Applying Mutation")
            # Mutation
            mutation_attempts, mutation_successes = self.apply_mutation(parents)

            #print("DEBUG: Applying Crossover")
            # Crossover
            crossover_attempts, crossover_successes = self.apply_crossover(parents)

            #print("DEBUG: Pruning Population")
            # Population Control
            self.prune_population()
            
            """
            Logging
            """
            #print("DEBUG: Logging Succes Rates")
            # Log sucess rates
            if mutation_attempts > 0:
                mutation_rate = (mutation_successes / mutation_attempts) * 100
                update_mutation_log(mutation_rate)
                update_latest_mutation_rates(mutation_rate)

            if crossover_attempts > 0:
                crossover_rate = (crossover_successes / crossover_attempts) * 100
                update_crossover_log(crossover_rate)
                update_latest_crossover_rates(crossover_rate)

            #print("DEBUG: Updating Global States")
            # Update global states
            population_size = sum(1 for node in self.working_population.nodes 
                                  if self.working_population.nodes[node]["level"] == "Individual")
            update_current_population_size(population_size)

            # TODO - this just doubles the space of working_population - I should probably
            #        just keep a global state of the top_n mols for the standings, rather than
            #        saving the entire population to global state. Or keep for reproduciblity?
            update_latest_population(self.working_population)

            #print("DEBUG: Updating Fitness Log")
            # Extract min, mean, and max fitness values
            individuals = [node for node in self.working_population.nodes 
                    if self.working_population.nodes[node]["level"] == "Individual"]
            fitness_values = [self.working_population.nodes[node]["raw_fitness"] 
                              for node in individuals]
            if fitness_values:
                min_fitness = min(fitness_values)
                mean_fitness = sum(fitness_values) / len(fitness_values)
                max_fitness = max(fitness_values)
                update_fitness_log((min_fitness, mean_fitness, max_fitness))

            #print("DEBUG: Updating Diversity Log")
            # Calculate and log diversity metrics
            diversity = calculate_population_diversity(self.working_population)
            self.diversity_log.append(diversity)
            self.working_population.nodes["Population"]["diversity"] = diversity
            update_latest_diversity(self.diversity_log)


        set_ga_active(False)
        return self.working_population, self.diversity_log
    
    def start_ga(self):
        # Load fragment library and pass to Fragment Mutation function
        # - These must be initialized here in order to do multi-processing batch runs
        # - otherwise they can cause serialization issues.
        

        """
        Starts execution of GA: initialize, run, log.
        """
        # Generate primary key run IDs
        run_id = str(uuid.uuid4())
        start_time = time.time()
        timestamp = datetime.datetime.now().isoformat()

        # Start tracking CPU usage
        process = psutil.Process()
        process.cpu_percent(interval=None)

        # Benchmarking and Logging
        #print("GA: GA complete, logging results")
        end_time = time.time()
        runtime = end_time - start_time

        # Initialize population and start the GA
        self.initialize_population()
        print("GA: Starting GA...")
        final_population, diversity_log = self.run_ga()

        # using oneshot here because I plan on running additional processing benchmarks
        with process.oneshot():
            cpu_usage = process.cpu_percent(interval=None)
            # TODO - implement continuous memory tracking

        # Collect fitness results TODO - append other fitness results
        fitness_results = {
            'fitness_log': get_fitness_history_log(),
            'mutation_rate_log': get_mutation_log(),
            'crossover_rate_log': get_crossover_log(),
            'diversity': self.diversity_log,
            }

        log_metrics(run_id,
                    timestamp,
                    runtime=runtime,
                    hyperparameters=self.config,
                    fitness_results=fitness_results,
                    validation_results={}, # TODO - handle val res
                    mem_usage=0, # TODO - handle mem calc
                    cpu_usage=cpu_usage,
                    population_type=self.population_type,
                    fitness_config = {"tuning_weights": get_tuning_weights(),
                                      "active_filters": get_active_filters(),
                                      "mw_range": get_mw_range(),
                                      "logp_range": get_logp_range()}
                    )
        print("GA: Simulation complete.")

if __name__ == "__main__":
    config = load_config()
    set_ga_active(True)
    ga = GeneticAlgorithm(config)
    final_population, diversity_log = ga.start_ga()