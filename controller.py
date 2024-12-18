from src.evolutionary_system.utils.config_loader import load_config
from src.evolutionary_system.utils.logger import log_metrics
from src.evolutionary_system.ga import run_ga
import datetime
import psutil



def start_ga():
    timestamp = datetime.datetime.now().isoformat()
    
    # TODO - generate unique run ID from timestamp
    
    config = load_config()
    
    # TODO - save/unpack current config - to be sent to logger.py
    
    stopped = False 

    # TODO - Start running metrics (e.g., runtime, memory, cpu usage)
    
    # TODO - call run_ga() w/ hyperparameters from config.yaml
    
    # TODO - handle end conditions (e.g, user defined or GA run completed)
    
    if stopped: 
        # TODO - Read in final fitness metrics from population
        # TODO - fill rest of parameters
        log_metrics(timestamp)