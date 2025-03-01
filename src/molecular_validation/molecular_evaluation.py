from guacamol.utils.chemistry import smiles_to_rdkit_mol #type:ignore
from guacamol.distribution_learning_benchmark import ValidityBenchmark, UniquenessBenchmark, NoveltyBenchmark, KLDivBenchmark
from moses.metrics.SA_Score.sascorer import calculateScore
import pandas as pd

class SimpleGenerator:
    """
    Wrapper class to make a list of SMILES compatible with Guacamol functions
    """
    def __init__(self, smiles_list):
        self.smiles_list = smiles_list
    def generate(self, number_samples):
        """
        Returns a subset of the SMILES list, mimicking a generator model.

        :param number_samples: Number of molecules to return.
        :return: List of SMILES strings.
        """
        return self.smiles_list[:number_samples]

# TODO add refereence mol set
def evaluate_mols(population, top_n=10, reference_set=None):
    """
    Evaluates top_n molecules from the population using GuacaMol.

    :param population: NetworkX Di.Graph of population.
    :param top_n: Number of top molecules to evaluate.
    :param reference_set: Set of reference molecules fo  calculating KL/FCD/Novelty of generated mols.
    :return: Pandas DataFrame with evaluation results.
    """
    try:
        
        # Extract top 10 Candidates from most recent population
        population_list = [
            (node, float(population.nodes[node].get("raw_fitness") or 0))
            for node in population.nodes
            if population.nodes[node]["level"] == "Individual"
        ]
        top_candidates = sorted(population_list, key=lambda x: x[1], reverse=True)[:top_n]

        results = []
        ga_output_set = []
        # Extract data from each Mol and append to results list, also calulating Synthetic Accessibility score
        for idx, (node, fitness) in enumerate(top_candidates):
            compounds = population.nodes[node]["compounds"]

            smiles = compounds[0]["structure"] if compounds else None

            if not smiles:
                continue

            mol = smiles_to_rdkit_mol(smiles)
            ga_output_set.append(smiles)

            if mol:
                sa_score = calculateScore(mol)

                results.append({
                    "Molecule": f"Compound {idx+1}",
                    "Fitness": round(fitness, 4),
                    "SAScore": round(sa_score, 4),
                    "SMILES": smiles
                })
        df = pd.DataFrame(results)

        # assess_model() expects a Generator model and exact number of samples stated
        ga_output_set = SimpleGenerator(ga_output_set)
        num_samples_available = len(ga_output_set.smiles_list)
        
        validity_benchmark = ValidityBenchmark(number_samples=num_samples_available)
        uniqueness_benchmark = UniquenessBenchmark(number_samples = num_samples_available)
        try:
            validity_score = validity_benchmark.assess_model(ga_output_set).score
            uniqueness_score = uniqueness_benchmark.assess_model(ga_output_set).score
        except Exception as e:
            print(f"Error : {e}")

        df["Validity"] = round(validity_score, 4)
        df["Uniqueness"] = round(uniqueness_score, 4)

        # Only calculate if reference set is present
        if reference_set:
            novelty_benchmark = NoveltyBenchmark(number_samples=num_samples_available, training_set=reference_set)
            kl_div_benchmark = KLDivBenchmark(number_samples=num_samples_available, training_set=reference_set)

            novelty_score = novelty_benchmark.assess_model(ga_output_set).score
            kl_div_score = kl_div_benchmark.assess_model(ga_output_set).score

            df["Novelty"] = round(novelty_score, 4)
            df["KL Divergence"] = round(kl_div_score, 4)
        return df


    
    except Exception as e:
        print("Error during GuacaMol Evaluation: {e}")
        return None