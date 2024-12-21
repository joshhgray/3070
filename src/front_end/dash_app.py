import dash #type: ignore
from dash import dcc, html, Input, Output, State #type: ignore
from dash.dependencies import Input, Output #type: ignore
import pandas as pd #type: ignore
import plotly.express as px #type: ignore
import yaml # type: ignore
import os
from threading import Thread
from controller import start_ga, stop_ga

app = dash.Dash(__name__)

config_path = os.path.abspath("config.yaml")

def load_hyperparameters():
    with open(config_path, "r") as file:
        return yaml.safe_load(file)

def save_hyperparameters(hyperparameters):
    with open(config_path, "w") as file:
        yaml.dump(hyperparameters, file)
    
config = load_hyperparameters()

app.layout = html.Div([
    html.H1("GA Dashboard"),

    html.Div([
        # Sliders for hyperparameters
        html.Label("Population Size"),
        dcc.Slider(id="population-size-slider", min=10, max=800, step=10, value=config.get("population_size")),
        html.Label("Number of Generations"),
        dcc.Slider(id="num-generations-slider", min=5, max=500, step=10, value=config.get("num_generations")),
        html.Label("Mutation Rate"),
        dcc.Slider(id="mutation-rate-slider", min=0.1, max=1.0, step=0.1, value=config.get("mutation_rate")),
        html.Label("Crossover Rate"),
        dcc.Slider(id="crossover-rate-slider", min=0.1, max=1.0, step=0.1, value=config.get("crossover_rate")),
        html.Label("Number of Elites (Individuals)"),
        dcc.Slider(id="num-elite-individuals-slider", min=1, max=10, step=1, value=config.get("num_elite_individuals")),
        html.Label("Number of Elites (Groups)"),
        dcc.Slider(id="num-elite-groups-slider", min=1, max=10, step=1, value=config.get("num_elite_groups")),
        html.Label("Number of threads"),
        # TODO - adjust when implementing Multi-threading
        dcc.Slider(id="num-threads-slider", min=1, max=1, step=1, value=1),

        # Dropdown for selection method
        html.Div([
            html.Label("Selection Method"),
            dcc.Dropdown(
                id="selection-method-dropdown",
                options=[
                    {"label": "Stochastic Universal Selection", "value": "sus"},
                    # TODO - Add other selection methods
                    ],
                value="sus",
                clearable=False
            )
        ]),

        # Inputs for fitness weights
        html.Div([
            html.Label("Fitness Weights"),
            html.Div([
                html.Label("Scaling Factor A"),
                dcc.Input(
                    id="fitness-weight-a",
                    type="number",
                    value=1.0, step=0.1
                )
            ]),
            html.Div([
                html.Label("Scaling Factor B"),
                dcc.Input(
                    id="fitness-weight-b",
                    type="number",
                    value=1.0, step=1
                )
            ]),
        ]),

        # Save Button
        html.Button("Save Config", id="save-button", n_clicks=0),
        html.Div(id="save-status")
    ]),

    html.Div([
        html.Button("Start GA", id="start-ga-button", n_clicks=0),
        html.Button("Stop GA", id="stop-ga-button", n_clicks=0),
        html.Div(id="ga-status")
    ]),

    # TODO - integrate graphics
    # placeholder
    dcc.Graph(id="ga-viz", style={"hegiht": "500px"}),
    
    dcc.Store(id="current-config", data=config),
    dcc.Store(id="ga-running", data=False),
    
])

@app.callback(
    Output("current-config", "data"),
    [
        Input("population-size-slider", "value"),
        Input("num-generations-slider", "value"),
        Input("mutation-rate-slider", "value"),
        Input("crossover-rate-slider", "value"),
        Input("num-elite-individuals-slider", "value"),
        Input("num-elite-groups-slider", "value"),
        Input("selection-method-dropdown", "value"),
        Input("fitness-weight-a", "value"),
        Input("fitness-weight-b", "value"),
        Input("num-threads-slider", "value")
    ],
    State("current-config", "data")
)
def update_config(population_size, 
                  num_generations, 
                  mutation_rate, 
                  crossover_rate, 
                  num_elite_individuals, 
                  num_elite_groups,
                  selection_method,
                  fitness_weight_a,
                  fitness_weight_b,
                  num_threads,
                  current_config):
    current_config["population_size"] = population_size
    current_config["num_generations"] = num_generations
    current_config["mutation_rate"] = mutation_rate
    current_config["crossover_rate"] = crossover_rate
    current_config["num_elite_individuals"] = num_elite_individuals
    current_config["num_elite_groups"] = num_elite_groups
    current_config["selection_method"] = selection_method
    current_config["fitness_weigths"] = [fitness_weight_a, fitness_weight_b]
    current_config["num_threads"] = num_threads
    return current_config

@app.callback(
    Output("save-status", "children"),
    Input("save-button", "n_clicks"),
    State("current-config", "data")
)
def save_config(n_clicks, current_config):
    if n_clicks > 0:
        save_hyperparameters(current_config)
        return "Saved config."
    return ""

@app.callback(
    [Output("ga-status", "children"),
     Output("ga-running", "data")],
     Input("start-ga-button", "n_clicks"),
     State("ga-running", "data")
)
def start_ga_callback(n_clicks, is_running):
    if n_clicks > 0:
        if not is_running:
            # open new thread to run the GA
            thread = Thread(target=start_ga)
            thread.start()
            return "GA Started..."
        return "GA already running."
    return "Click 'Start GA'.", is_running

@app.callback(
    [Output("ga-status", "children"),
     Output("ga-running", "data")],
     Input("stop-ga-button", "n_clicks"),
     State("ga-running", "data")
)
def stop_ga_callback(n_clicks, is_running):
    if n_clicks > 0:
        if is_running:
            stop_ga()
            return "Stopping GA..."
        return "GA not running."
    return "Click 'Stop GA' to halt the process.", is_running



if __name__ == "__main__":
    app.run_server(debug=True)

