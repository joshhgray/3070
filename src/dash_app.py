import dash
from dash import dcc, html, Input, Output, State, dash_table
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
from src.molecular_validation.molecular_evaluation import evaluate_mols
from src.evolutionary_system.utils.ga_state import get_latest_population
from src.evolutionary_system.utils.nx_graph_to_mol import nx_graph_to_mol
from threading import Thread
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import base64
from src.data_pipeline.mol_to_graph import mol_to_graph


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.DARKLY], suppress_callback_exceptions=True)

config_path = os.path.abspath("config.yaml")
last_diversity_fig = None

# TODO - move functions to utils? idk
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
                dcc.Input(id="pop-size", type="number", value=1000, min=10, max=2000,
                        style={
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
                dcc.Input(id="num-gen", type="number", value=1000, min=1, max=9000,
                        style={
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
                dbc.Label("Carrying Capacity"),
                dcc.Input(id="carrying-capacity", type="number", value=2000, min=10, max=5000,
                        style={
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
            dbc.Label("Selection Methods"),
            dcc.Dropdown(
                id="selection-method",
                options=[
                    {"label": "SUS", "value": "stochastic_universal_sampling"},
                    {"label": "Ranked", "value": "rank_based_selection"}
                ],
                value="stochastic_universal_sampling"
            )
        ], className="mb-3"),
        html.Div([
            dbc.Label("Mutation Methods"),
            dcc.Dropdown(
                id="mutation-methods",
                options=[
                    {"label": "Hydroxylate", "value": "hydroxylate_mutate"},
                    {"label": "Atomic Substitution", "value": "atomic_substitution"}
                ],
                value=["hydroxylate_mutate", "atomic_substitution"],
                multi=True
            )
        ], className="mb-3"),
        html.Div([
            dbc.Label("Crossover Methods"),
            dcc.Dropdown(
                id="crossover-methods",
                options=[
                    {"label": "Graph Based", "value": "graph_based_crossover"},
                    #{"label": "HGT", "value": "hybrid_gene_crossover"}
                ],
                value=["graph_based_crossover"], 
                multi=True 
            )
        ], className="mb-3"),
        html.Div([
            dbc.Label("Fitness Function"),
            dcc.Dropdown(
                id="fitness-functions",
                options=[
                    {"label": "QED", "value": "calculate_qed"},
                    #{"label": "Ro5", "value": "rule_of_five"}
                ],
                value=["calculate_qed"],
                multi=True
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
live_tracking_panel = dbc.Card([
    dbc.CardHeader("Live Tracking"),
    dbc.CardBody([
        dbc.Row([
            dbc.Col(dcc.Graph(id="diversity-graph", style={"height": "375px", "width": "100%", "height": "375px", "padding": "0", "margin": "0"}), width=3),
            dbc.Col(dcc.Graph(id="graph-2", style={"height": "375px", "width": "100%", "height": "375px", "padding": "0", "margin": "0"}), width=3),
            dbc.Col(dcc.Graph(id="graph-3", style={"height": "375px", "width": "100%", "height": "375px", "padding": "0", "margin": "0"}), width=3),
            dbc.Col(dcc.Graph(id="graph-4", style={"height": "375px", "width": "100%", "height": "375px", "padding": "0", "margin": "0"}), width=3)
        ])
    ]),
    ],
    className="mb-3"
)

# Results Panel
results_panel = dbc.Card([
    dbc.CardHeader("Top 10 Compounds"),
    dbc.CardBody(id="results-panel-body", children=[
        html.P("Standings will be updated here.")
    ])
], className="mb-3")

molecule_image_section = dbc.Card([
    dbc.CardHeader("Molecular Structure"),
    dbc.CardBody([
        html.Img(id="molecule-image", style={"width": "300px", "height": "300px"}),
        html.P(id="molecule-name")
    ])
], className="mb-3")

app.layout = dbc.Container([
    dcc.Store(id="molecule-list"),
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
        dbc.Col([
            dbc.Row(dbc.Col(live_tracking_panel), className="mb-3"),
            dbc.Row([
                dbc.Col(results_panel, width=7),
                dbc.Col(molecule_image_section, width=3)
            ])
        ], width=10)
    ])
], fluid=True)

@app.callback(
    Output("current-config", "data"),
    [
        Input("pop-size", "value"),
        Input("num-gen", "value"),
        Input("carrying-capacity", "value"),
        Input("selection-method", "value"),
        Input("mutation-methods", "value"),
        Input("crossover-methods", "value"),
        Input("fitness-functions", "value"),
    ],
    State("current-config", "data")
)
def update_config(population_size, num_generations, carrying_capacity, selection_method,
                  mutation_methods, crossover_methods, fitness_functions, current_config):
    # Slider updates result in the config converting to a string
    # So, it must be manually converted back to dict here
    current_config = {}
    current_config["population_size"] = population_size
    current_config["num_generations"] = num_generations
    current_config["carrying_capacity"] = carrying_capacity
    current_config["selection_method"] = selection_method
    current_config["mutation_methods"] = mutation_methods
    current_config["crossover_methods"] = crossover_methods
    current_config["fitness_functions"] = fitness_functions

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
    # TODO - figure out how to hide default graph from 
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
            labels={"x": "", "y": ""}
    )

    # Remove margins and disable tick labels
    last_diversity_fig.update_layout(
        autosize=True,
        margin=dict(l=0,r=0,t=0,b=0,pad=0),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)"
    )

    last_diversity_fig.update_layout(
        annotations=[
            dict(
                text="Population Diversity",
                x=0.5,
                y=0.5,
                xref="paper",
                yref="paper",
                font=dict(size=30, color="rgba(255,255,255,0.2)"),
                showarrow=False
            ),
            
        ],
        #autosize=True,
    )
    last_diversity_fig.update_xaxes(
        showticklabels=False,
        ticks="",
        fixedrange=True,  # Prevents zooming
        showgrid=False,
        zeroline=False,
    )
    last_diversity_fig.update_yaxes(
        showticklabels=True,
        ticks="outside",
        tickvals=[0.0,0.5,1.0],
        tickfont=dict(size=12),
        range=[0,1],
        fixedrange=True,
        showgrid=False,
        zeroline=False,
    )
    return last_diversity_fig

@app.callback(
    Output("results-panel-body", "children"),
    Output("molecule-list", "data"),
    Input("evaluate-btn", "n_clicks"),
    prevent_initial_call=True
)
def update_standings(n_clicks):
    if n_clicks:
        try:
            population = get_latest_population()
            # TODO - find a reference set (CheMBL)
            reference_set = None
            df = evaluate_mols(population, top_n=10, reference_set=reference_set)
            mol_list = df.to_dict("records")

            if df is not None:
                return [
                    dash_table.DataTable(
                        id="evaluation-table",
                        columns=[{"name": col, "id": col} for col in df.columns],
                        data=mol_list,
                        style_table={"overflowX": "auto"},
                        style_cell={"textAlign": "center", "color": "white", "backgroundColor": "#2a2a2a", "padding": "5px"},
                        style_header={"backgroundColor": "#1f1f1f", "fontWeight": "bold", "color": "white"},
                        sort_action="native",
                        row_selectable="single",
                        page_size=10
                    )
                ], mol_list
            else:
                return [html.P("Error evaluating molecule.")], None
            
        except Exception as e:
            print(f"Error Loading Standings: {e}")
            return [html.P(f"Error Loading Standings")], None
    return [html.P("Click 'Evaluate Current Population' to view standings")], None

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

@app.callback(
    Output("molecule-image", "src"),
    Output("molecule-name", "children"),
    Input("evaluation-table", "selected_rows"),
    State("molecule-list", "data"),
    prevent_initial_call=True
)
def update_molecule_image(selected_rows, mol_list):
    if not selected_rows or mol_list is None:
        return dash.no_update, dash.no_update

    try:
        selected_index = selected_rows[0]
        selected_mol = mol_list[selected_index]
        smiles = selected_mol.get("SMILES")

        if not smiles:
            return dash.no_update, selected_mol.get("Molecule", "Error Loading Compound Image")

        mol_graph = mol_to_graph(smiles)
        mol = nx_graph_to_mol(mol_graph)

        if not mol:
            return dash.no_update, selected_mol.get("Molecule", "Error Loading Compound Image")
        
        # open buffer and store image
        img = Draw.MolToImage(mol, size=(300, 300))
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        # Bytes -> B64 -> UTF-8
        encoded_image = base64.b64encode(buffered.getvalue()).decode("utf-8")

        return f"data:image/png;base64, {encoded_image}", selected_mol.get("Molecule", "Error Loading Compound Image")
        
    except Exception as e:
        print(f"Error updating Molecular Image: {e}")
        return dash.no_update, "Error displaying molecular image."
    
if __name__ == "__main__":
    app.run_server(debug=True)