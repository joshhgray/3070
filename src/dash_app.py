import dash #type: ignore
from dash import dcc, html, Input, Output, State #type: ignore
from dash.dependencies import Input, Output#type: ignore
import dash_bootstrap_components as dbc #type: ignore
import pandas as pd #type: ignore
import plotly.express as px #type: ignore
import yaml # type: ignore
import os
from threading import Thread
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from .controller import start_ga#, stop_ga
import numpy as np # type: ignore

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.YETI])

config_path = os.path.abspath("config.yaml")

def load_hyperparameters():
    with open(config_path, "r") as file:
        return yaml.safe_load(file)

def save_hyperparameters(hyperparameters):
    with open(config_path, "w") as file:
        yaml.dump(hyperparameters, file)

def load_diversity_log(file_path="src/evolutionary_system/utils/metrics_log.csv"):
    try:
        data = pd.read_csv(file_path)
        fitness_results = data["fitness_results"].iloc[-1]
        fitness_results = eval(fitness_results) 
        diversity_log = fitness_results.get("diversity", [])
        return [float(value) for value in diversity_log]
    except (FileNotFoundError, KeyError, SyntaxError) as e:
        print(f"Error loading diversity log: {e}")
        return []
    
config = load_hyperparameters()

app.layout = dbc.Container([
    dbc.Row(
        dbc.Col(html.H1("GA Dashboard"), width={"size": 4, "offset": 3}),
        className="mb-4"
    ),
    dbc.Row([
        # Config
        dbc.Col([
            html.H2("Configuration"),
                # Hyperparameters Adjustment Section
                dbc.Card([
                    dbc.CardHeader("Hyperparameters"),
                    dbc.CardBody([
                        # Population Size Slider
                        dbc.Row([
                            dbc.Col(dbc.Label("Population Size", html_for="population-size-slider"), width=2),
                            dbc.Col(dcc.Slider(id="population-size-slider", min=50, max=800, step=50, value=config.get("population_size")), width=6),
                        ]),
                        # Number of Generations Slider
                        dbc.Row([
                            dbc.Col(dbc.Label("Number of Generations", html_for="num-generations-slider"), width=2),
                            dbc.Col(dcc.Slider(id="num-generations-slider", min=5, max=500, step=50, value=config.get("num_generations")), width=6),
                        ]),
                        # Mutation Rate Slider
                        dbc.Row([
                            dbc.Col(dbc.Label("Mutation Rate", html_for="mutation-rate-slider"), width=2),
                            dbc.Col(dcc.Slider(id="mutation-rate-slider", min=0.1, max=1.0, step=0.1, value=config.get("mutation_rate")), width=6),
                        ]),
                        # Crossover Rate Slider
                        dbc.Row([
                            dbc.Col(dbc.Label("Crossover Rate", html_for="crossover-rate-slider"), width=2),
                            dbc.Col(dcc.Slider(id="crossover-rate-slider", min=0.1, max=1.0, step=0.1, value=config.get("crossover_rate")), width=6),
                        ], className="conf"),
                        # Number of Elite Individuals Slider
                        dbc.Row([
                            dbc.Col(dbc.Label("Number of Elites (Individuals)", html_for="num-elite-individuals-slider"), width=2),
                            dbc.Col(dcc.Slider(id="num-elite-individuals-slider", min=1, max=10, step=1, value=config.get("num_elite_individuals")), width=6),
                        ]),
                        # Number of Elite Groups Slider
                        dbc.Row([
                            dbc.Col(dbc.Label("Number of Elites (Groups)", html_for="num-elite-groups-slider"), width=2),
                            dbc.Col(dcc.Slider(id="num-elite-groups-slider", min=1, max=5, step=1, value=config.get("num_elite_groups")), width=6),
                        ]),
                        # Number of Threads Slider
                        # TODO - adjust when implementing Multi-threading
                        dbc.Row([
                            dbc.Col(dbc.Label("Number of Threads", html_for="num-threads-slider"), width=2),
                            dbc.Col(dcc.Slider(id="num-threads-slider", min=1, max=1, step=1, value=1), width=6),
                        ]),
                        # Selection Method Dropdown Menu
                        dbc.Row([
                            dbc.Col(dbc.Label("Selection Method", html_for='selection-method-dropdown'), width=2),
                            dbc.Col(dcc.Dropdown(id='selection-method-dropdown', options=[{"label": "Stochastic Universal Sampling", "value": "sus"}], value="sus", clearable=False))
                        ]),
                        # Fitness Weight Selection
                        dbc.Row([
                            dbc.Col(dbc.Label("Fitness Weights")),
                            dbc.Col(dbc.Label("Scaling Factor A", html_for='fitness-weight-a'), width=1),
                            dbc.Col(dcc.Input(id="fitness-weight-a", type="number", value=1.0, step=0.1)),
                            dbc.Col(dbc.Label("Scaling Factor B", html_for='fitness-weight-b'), width=1),
                            dbc.Col(dcc.Input(id="fitness-weight-b", type="number", value=1.0, step=0.1)),

                        ]),
                        
                    ]),
                ]),
                dbc.Button("Start GA", id="start-ga-button", color="primary"),
                dbc.Button("Stop GA", id="stop-ga-button", color="primary"),
                dcc.Store(id="current-config", data=config),
                dcc.Store(id="ga-running", data=False),
                dbc.Badge(id="ga-status"), 

                # Visualizations Sections
                dbc.Col([                     
                    html.H2("Visualizations"),
                    dbc.Button("Update Results", id="update-results-button"),
                    dbc.Card([
                        dbc.CardHeader("Population Level Diversity by Generation"),
                        dbc.CardBody(dcc.Graph(id="diversity-graph")),
                    ]),
                ], width=8),
        ]),
])
], fluid=True),


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
    current_config["fitness_weights"] = [fitness_weight_a, fitness_weight_b]
    current_config["num_threads"] = num_threads

    save_hyperparameters(current_config)
    return "Configuration Updated."


@app.callback(
    [Output("ga-status", "children", allow_duplicate=True),
     Output("ga-running", "data", allow_duplicate=True)],
     Input("start-ga-button", "n_clicks"),
     State("ga-running", "data"),
     # Initial call prevention required for duplicate outputs
     prevent_initial_call=True
)
def start_ga_callback(n_clicks, is_running):
    if n_clicks and n_clicks > 0:
        if not is_running:
            # open new thread to run the GA
            thread = Thread(target=start_ga)
            thread.start()
            return "GA Started...", is_running
        return "GA already running.", is_running
    return "Click 'Start GA'.", is_running

# TODO - Uncomment once fully implemented in backend
@app.callback(
    [Output("ga-status", "children", allow_duplicate=True),
     Output("ga-running", "data", allow_duplicate=True)],
     Input("stop-ga-button", "n_clicks"),
     State("ga-running", "data"),
     # Initial call prevention required for duplicate outputs
     prevent_initial_call=True,
)
def stop_ga_callback(n_clicks, is_running):
    if n_clicks and n_clicks > 0:
        return "", is_running
    return "", is_running


@app.callback(
        Output("diversity-graph", "figure"),
        Input("update-results-button", "n_clicks")
)
def update_diversity_graph(n_clicks):
    if n_clicks and n_clicks > 0:
        diversity_log = load_diversity_log()
        generations = list(range(1, len(diversity_log) + 1))
        fig = px.line(
            x=generations,
            y=diversity_log,
            labels={"x": "Generation", "y": "Diversity"},
            title="Diversity by Generation"
        )
        return fig
    return px.line(
        labels={"x": "Generation", "y": "Diversity"},
        title="Diversity by Generation"
    )


if __name__ == "__main__":
    app.run_server(debug=True)