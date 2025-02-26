import dash
from dash import dcc, html, Input, Output, State
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
import yaml
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from .controller import start_ga
from src.evolutionary_system.utils.ga_state import set_ga_active, is_ga_active, get_latest_diversity, get_latest_population
from threading import Thread

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.DARKLY])

config_path = os.path.abspath("config.yaml")
last_diversity_fig = None

# TODO - move functions to utils?
def load_hyperparameters():
    with open(config_path, "r") as file:
        return yaml.safe_load(file)

def save_hyperparameters(hyperparameters):
    with open(config_path, "w") as file:
        yaml.dump(hyperparameters, file)

def load_diversity_log(file_path="src/evolutionary_system/utils/metrics_log.csv"):
    try:
        data = pd.read_csv(file_path)
        return list(data["diversity_score"])
    except Exception as e:
        print(f"Error loading diversity log: {e}")
        return []

config = load_hyperparameters()

# Header
header = dbc.Row(
    dbc.Col(html.H1("Genetic Algorithm Dashboard"), width={"size": 6}),
    className="mb-4"
)

# General Configuration Panel
general_config_panel = dbc.Card([
    dbc.CardHeader("General Configuration"),
    dbc.CardBody([
        dbc.Row([
            dbc.Col([
                dbc.Label("Population Size"),
                # Max population set to 2000 currently due to size of data present
                dcc.Input(id="pop-size", type="number", value=100,
                    min=10, max=2000, style={
                        "width": "110px",
                        "textAlign": "center",
                        "borderRadius": "8px",
                        "border": "1px solid #ccc",
                        "appearance": "textfield",
                        "padding": "5px"
                    }
                )
            ]),
            dbc.Col([
                dbc.Label("Generations"),
                # Artificial cap of 9000 generations - this should be bypassable in advanced settings
                dcc.Input(id="num-gen", type="number", value=50,
                    min=1, max=9000, style={
                        "width": "110px",
                        "textAlign": "center",
                        "borderRadius": "8px",
                        "border": "1px solid #ccc",
                        "appearance": "textfield",
                        "padding": "5px"
                    }
                )
            ])
        ])
    ])
], className="mb-3")

# Genetic Algorithm Operators Panel
ga_operators_panel = dbc.Card([
    dbc.CardHeader("Genetic Algorithm Operators"),
    dbc.CardBody([
        html.Div([
            dbc.Label("Selection Method"),
            dcc.Dropdown(
                id="selection-method",
                options=[
                    {"label": "Verhulst", "value": "verhulst"},
                    {"label": "Ranked", "value": "ranked"}
                ],
                value="ranked"
            )
        ], className="mb-3"),
        html.Div([
            dbc.Label("Mutation Method"),
            dcc.Dropdown(
                id="mutation-method",
                options=[
                    {"label": "Hydroxylate", "value": "hydroxylate"},
                    {"label": "Bit Flip", "value": "bitflip"}
                ],
                value="hydroxylate"
            )
        ], className="mb-3"),
        html.Div([
            dbc.Label("Crossover Method"),
            dcc.Dropdown(
                id="crossover-method",
                options=[
                    {"label": "One-Point", "value": "onepoint"},
                    {"label": "Two-Point", "value": "twopoint"}
                ],
                value="onepoint"
            )
        ], className="mb-3"),
        html.Div([
            dbc.Label("Fitness Function"),
            dcc.Dropdown(
                id="fitness-function",
                options=[
                    {"label": "QED", "value": "qed"},
                    {"label": "Ro5", "value": "rule_of_five"}
                ],
                value="qed"
            )
        ], className="mb-3")
    ])
], className="mb-3")

# Control Buttons Panel
control_buttons_panel = dbc.Card([
    dbc.CardHeader("GA Control"),
    dbc.CardBody([
        dbc.Button("Start GA", id="start-btn", color="success", className="mb-2", size="lg"),
        dbc.Button("Stop GA", id="stop-btn", color="danger", className="mb-2", size="lg"),
        dbc.Button("Evaluate Current Population", id="evaluate-btn", color="primary", className="mb-2", size="md")
    ])
], className="mb-3")

# Live Tracking Panel (right side)
live_tracking_panel = dbc.Card(
    dbc.CardBody([
        html.H3("Live Tracking"),
        dcc.Graph(id="diversity-graph")
    ]),
    className="mb-3"
)

# Results Panel
results_panel = dbc.Card([
    dbc.CardHeader("Top 10 Compounds"),
    dbc.CardBody(id="results-panel-body", children=[
        html.P("Standings will be updated here.")
    ])
], className="mb-3")

app.layout = dbc.Container([
    # Store current config as default
    dcc.Store(id="current-config", data=load_hyperparameters()),
    dcc.Store(id="ga-running", data=False),
    # Interval to update live tracking every 1 second also used for checking ga status
    dcc.Interval(
        id="interval-component",
        interval=1000,
        n_intervals=0
    ),
    header,
    dbc.Row([
        # Left column for configuration panels
        dbc.Col([
            general_config_panel,
            ga_operators_panel,
            control_buttons_panel
        ], width=2),
        # Middle column, live tracking
        dbc.Col(live_tracking_panel, width=7),
        # right side column, results and standings
        dbc.Col(results_panel, width=3)
    ])
], fluid=True)

@app.callback(
    Output("current-config", "data"),
    [
        Input("pop-size", "value"),
        Input("num-gen", "value"),
        Input("selection-method", "value"),
        Input("mutation-method", "value"),
        Input("crossover-method", "value"),
        Input("fitness-function", "value"),
    ],
    State("current-config", "data")
)
def update_config(population_size, num_generations, selection_method,
                  mutation_method, crossover_method, fitness_function,
                  current_config):
    # Slider updates result in the config converting to a string
    # So, it must be manually converted back to dict here
    current_config = {}
    current_config["population_size"] = population_size
    current_config["num_generations"] = num_generations
    current_config["selection_method"] = selection_method
    current_config["mutation_method"] = mutation_method
    current_config["crossover_method"] = crossover_method
    current_config["fitness_function"] = fitness_function



    save_hyperparameters(current_config)
    return "Configuration Updated."

@app.callback(
    Output("ga-running", "data", allow_duplicate=True),
     Input("start-btn", "n_clicks"),
     State("ga-running", "data"),
     # Initial call prevention required for duplicate outputs 
     prevent_initial_call=True
)
def start_ga_callback(n_clicks, is_running):
    if n_clicks:
        if not is_running:
            set_ga_active(True)
            # open new thread to run the GA
            thread = Thread(target=start_ga)
            thread.start()
            return True
        return True
    return False

@app.callback(
    Output("ga-running", "data", allow_duplicate=True),
     Input("stop-btn", "n_clicks"),
     State("ga-running", "data"),
     # Initial call prevention required for duplicate outputs
     prevent_initial_call=True,
)
def stop_ga_callback(n_clicks, is_running):
    if n_clicks and is_running:
        set_ga_active(False)
        return False
    return is_running

@app.callback(
        Output("diversity-graph", "figure", allow_duplicate=True),
        Input("interval-component", "n_intervals"),
        State("ga-running", "data"),
        # Required for duplicate outputs
        prevent_initial_call=True
)
def update_diversity_graph(n_intervals, is_running):
    global last_diversity_fig

    # Displays empty graph when GA hasn't been run, or displays last frame when GA is stopped
    if not is_ga_active():
        return last_diversity_fig if last_diversity_fig else px.line(title="Population Diversity")
    
    diversity_log = get_latest_diversity()

    # Prevent update if there is no log yet, or display most recent figure
    if not diversity_log:
        print("dash_app: No diversity score log found.")
        return last_diversity_fig if last_diversity_fig else dash.no_update
    
    scores = diversity_log
    generations = list(range(1, len(diversity_log) + 1))
    last_diversity_fig = px.line(
            x=generations,
            y=scores,
            labels={"x": "Generation", "y": "Diversity"},
            title="Diversity by Generation"
    )
    return last_diversity_fig

@app.callback(
    Output("results-panel-body", "children"),
    Input("evaluate-btn", "n_clicks"),
    prevent_initial_call=True
)
def update_standings(n_clicks):
    if n_clicks:
        try:
            population = get_latest_population()

            # Extract relevant information from population graph to list
            population_list = [
                (node, float(population.nodes[node].get("raw_fitness") or 0))
                for node in population.nodes
                if population.nodes[node]["level"] == "Individual"
            ]

            # Sort by fitness descending, extract top 10
            top_10_candidates = sorted(population_list, key=lambda x: x[1], reverse=True)[:10]

            standings = [
                html.P(f"{idx+1}. Fitness: {fitness:.4f}")
                for idx, (node, fitness) in enumerate(top_10_candidates)
            ]

            return standings

        except Exception as e:
            print(f"Error Loading Standings: {e}")
            return [html.P(f"Error Loading Standings")]
    return [html.P("Click 'Evaluate Current Population' to view standings")]

@app.callback(
    Output("ga-running", "data"),
    Input("interval-component", "n_intervals"),
    Input("ga-running", "data")
)
def update_ga_status(n_intervals, current_status):
    """
    Continuous polling of GA status to ensure synchronization with global state of GA 
    
    :param n_intervals: 1 second interval, this function will run every 1 second.
    :param current_status: Local status 
    """
    ga_status = is_ga_active()

    if ga_status != current_status:
        return ga_status
    return current_status

if __name__ == "__main__":
    app.run_server(debug=True)