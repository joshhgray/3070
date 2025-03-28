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
from src.evolutionary_system.GeneticAlgorithm import GeneticAlgorithm
from src.evolutionary_system.utils.config_loader import load_config
from src.evolutionary_system.utils.ga_state import (
    set_ga_active, is_ga_active, get_latest_diversity, get_latest_population, 
    get_latest_crossover_rates, get_latest_mutation_rates, get_crossover_log, 
    get_mutation_log, get_current_generation_number, get_current_population_size,
    get_selected_fitness_functions, update_selected_fitness_functions, update_mutation_probabilities, 
    get_fitness_history_log, update_active_filters, update_mw_range, update_logp_range, update_tuning_weights)
from src.molecular_validation.molecular_evaluation import evaluate_mols
from src.evolutionary_system.utils.ga_state import get_latest_population
from threading import Thread
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import base64


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
                # Max population artificial cap at 100k
                dcc.Input(id="pop-size", type="number", value=1000, min=10, max=100000,
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
                # Artificial cap of 24000 generations
                dcc.Input(id="num-gen", type="number", value=1000, min=1, max=24000,
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
                dcc.Input(id="carrying-capacity", type="number", value=2000, min=10, max=100000,
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

# Dataset Toggle Selection Panel
dataset_toggle_panel = dbc.Card([
    dbc.CardHeader("Data Source"),
    dbc.CardBody([
        dbc.Row([
            dbc.RadioItems(
                id="data-toggle",
                options=[
                    {"label": "Precompiled BGC Dataset", "value": "bgc"},
                    {"label": "Custom SMILES Dataset", "value": "smiles"},
                ],
                value='bgc'
            )
        ])
    ])
], className="mb-3")

# Fitness Function Panel
fitness_panel = dbc.Card([
    dbc.CardHeader("Fitness Function"),
    dbc.CardBody([
        html.P("Tuning: mean of selected.", style={"color": "grey"}),
        html.P("Filters: hard exclusion.", style={"color": "grey"}),
        dbc.Checklist(
            options=[{"label": "Advanced Mode", "value": "advanced"}],
            value=[],
            id="fitness-advanced-toggle",
            inline=True,
            style={"marginBottom": "10px"}
        ),
        dbc.Collapse(
            dbc.Checklist(
                id="fitness-function-selection",
                options=[
                    {"label": "QED", "value": "qed"},
                    {"label": "SA score", "value": "sa"},
                    #{"label": "Lipinski Score (Ro5)", "value": "ro5"},
                ],
                value=get_selected_fitness_functions(),
                inline=True,
            ),
            id="basic-fitness-collapse",
            is_open=True
        ),
        dbc.Collapse([
            # TUNING WEIGHTS
            html.H6("Tuning Weights", className="mt-3"),

            # QED
            html.Label("QED"),
            dcc.Slider(id="qed-tuning", min=0, max=1, value=0.5),

            # Synthetic Accessibility
            html.Label("Synthetic Accessibility (SA)"),
            dcc.Slider(id="sa-tuning", min=0, max=1, value=0.5),

            # Molecular Weight
            html.Label("Molecular Weight"),
            dcc.Slider(id="mw-tuning", min=0, max=1, value=0.5),
            html.Br(),

            # FILTER SETTINGS
            html.H6("Filter Settings"),

            # Molecular Weight Filter
            html.Label("Molecular Weight Range Filter"),
            dbc.InputGroup([
                dbc.InputGroupText("Min"),
                dbc.Input(id="mw-filter-min", type="number", value=250),
                dbc.InputGroupText("Max"),
                dbc.Input(id="mw-filter-max", type="number", value=500)
            ]),

            # LogP Filter
            html.Label("LogP Range Filter"),
            dbc.InputGroup([
                dbc.InputGroupText("Min"),
                dbc.Input(id="logp-filter-min", type="number", value=-2),
                dbc.InputGroupText("Max"),
                dbc.Input(id="logp-filter-max", type="number", value=7)
            ]),

            # Lipinski Ro5 Filter
            dbc.Checklist(
                id="filter-selection",
                options=[
                    {"label": "Lipinski Rule of Five", "value": "ro5"},
                ],
                value=["ro5"]
            )

        ], id="advanced-fitness-collapse", is_open=False),

        # ELITISM WEIGHT
        html.Br(),
        html.Div([
            html.Label("Elitism Weight"),
            dcc.Slider(id="elitism-weight", min=0, max=5, step=1, value=1,
                       marks={0: "0%", 1: "1%", 2: "2%", 3: "3%", 4: "4%", 5: "5%"}),
        ]),
        
    ])
], className="mb-3")

# Selection Choice
selection_panel = dbc.Card([
    dbc.CardHeader("Selection Method"),
    dbc.CardBody([
        html.Div([
            dcc.Dropdown(
                id="selection-method",
                options=[
                    {"label": "SUS", "value": "stochastic_universal_sampling"},
                    {"label": "Ranked", "value": "rank_based_selection"}
                ],
                value="stochastic_universal_sampling"
            )
        ], className="mb-3"),
    ])
], className="mb-3")

# Mutation Choices and Weight Sliders
mutation_panel = dbc.Card([
    dcc.Store(id="mutation-weights-store"),
    dbc.CardHeader("Mutation Type Weights"),
    dbc.CardBody([
        # MUTATION METHOD SELECTION
        html.Div([
            dbc.CardBody([
                html.Label("Atomic Substitution"),
                dcc.Slider(id="atomic-input", min=0, max=1, step=0.01, value=0.25,
                           marks={0: "0", 0.25: "0.25", 0.5: "0.5", 0.75: "0.75", 1: "1"}),

                html.Label("Functional Group Mutation"),
                dcc.Slider(id="functional-group-input", min=0, max=1, step=0.01, value=0.25,
                           marks={0: "0", 0.25: "0.25", 0.5: "0.5", 0.75: "0.75", 1: "1"}),

                html.Label("Ring Mutation"),
                dcc.Slider(id="ring-input", min=0, max=1, step=0.01, value=0.0,
                           marks={0: "0", 0.25: "0.25", 0.5: "0.5", 0.75: "0.75", 1: "1"}),

                html.Label("Fragment Mutation"),
                dcc.Slider(id="fragment-input", min=0, max=1, step=0.01, value=0.0,
                           marks={0: "0", 0.25: "0.25", 0.5: "0.5", 0.75: "0.75", 1: "1"}),
            ])
            ]),
        ]),
], className="mb-3")

# Crossover Choices
crossover_panel = dbc.Card([
    dbc.CardHeader("Crossover Method(s)"),
    dbc.CardBody([
        # CROSSOVER METHOD SELECTION
        html.Div([
            dcc.Dropdown(
                id="crossover-methods",
                options=[
                    {"label": "MCS Crossover", "value": "mcs_crossover"},
                    #{"label": "Graph Based", "value": "graph_based_crossover"},
                    #{"label": "HGT", "value": "hybrid_gene_crossover"}
                ],
                value=["mcs_crossover"], 
                multi=True 
            )
        ]),
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
            # LIVE TRACKING SLOT 1 - General GA Stats
            dbc.Col(
                html.Div([
                dbc.Row([
                    # Live Generation Number Display
                    dbc.Col([
                        dbc.Card([
                            dbc.CardBody([
                                html.H4(id="current-generation-number-display", className="card-title", style={"textAlign": "center"}),
                                html.P("Current Generation", className="card-text", style={"textAlign": "center"})
                            ])
                        ], style={"marginBottom": "10px", "backgroundColor": "#99b3ff"})
                    ]),
                    # Live Population Size Dislpay
                    dbc.Col([
                        dbc.Card([
                            dbc.CardBody([
                                html.H4(id="current-population-size-display", className="card-title", style={"textAlign": "center"}),
                                html.P("Current Population Size", className="card-text", style={"textAlign": "center"})
                            ])
                        ], style={"marginBottom": "10px", "backgroundColor": "#99b3ff"})
                    ])
                ]),
                dbc.Row([
                    # Live Mutation Success Rate Display
                    dbc.Col([
                        dbc.Card([
                            dbc.CardBody([
                                html.H5(id="mutation-percentage-display", className="card-title", style={"textAlign": "center"}),
                                html.P("Current Mutation Success Rate", className="card-text", style={"textAlign": "center"})
                            ])
                        ], style={"marginBottom": "10px", "backgroundColor": "#809fff"})
                    ]),
                    # Live Average Mutation Success Rate Display
                    dbc.Col([
                        dbc.Card([
                            dbc.CardBody([
                                html.H5(id="avg-mutation-percentage-display", className="card-title", style={"textAlign": "center"}),
                                html.P("Average Mutation Success Rate", className="card-text", style={"textAlign": "center"})
                            ])
                        ],style={"marginBottom": "10px", "backgroundColor": "#809fff"})
                    ])
                ]),
                dbc.Row([
                    # Live Crossover Success Rate Display
                    dbc.Col([
                        dbc.Card([
                            dbc.CardBody([
                                html.H5(id="crossover-percentage-display", className="card-title", style={"textAlign": "center"}),
                                html.P("Current Crossover Success Rate", className="card-text", style={"textAlign": "center"})
                            ])
                        ], style={"marginBottom": "10px", "backgroundColor": "#668cff"})
                    ]),
                    # Live Average Crossover Success Rate Display
                    dbc.Col([
                        dbc.Card([
                            dbc.CardBody([
                                html.H5(id="avg-crossover-percentage-display", className="card-title", style={"textAlign": "center"}),
                                html.P("Average Crossover Success Rate", className="card-text", style={"textAlign": "center"})
                            ])
                        ], style={"marginBottom": "10px", "backgroundColor": "#668cff"})
                    ])
                ]),
            ]), width=4),
            # LIVE TRACKING SLOT 2 - Population Diversity Graph
            dbc.Col(dcc.Graph(id="diversity-graph", style={"height": "100%", "width": "100%", "padding": "0", "margin": "0"}), width=4),
            # LIVE TRACKING SLOT 3
            dbc.Col(dcc.Graph(id="graph-3", style={"height": "100%", "width": "100%", "padding": "0", "margin": "0"}), width=4),
        ])
    ]),
],className="mb-3")

# Results Panel
results_panel = dbc.Card([
    dbc.CardHeader("Top 10 Compounds"),
    dbc.CardBody(
        id="results-panel-body", 
        children=[html.P("Standings will be updated here.")]
    )
], className="mb-3")

molecule_image_section = dbc.Card([
    dbc.CardHeader("Molecular Structure"),
    dbc.CardBody([
        html.Div(
            html.Img(id="molecule-image", style={"width": "300px", "height": "300px"}),
            id="molecule-image-container",
            # centers the image
            style={"textAlign": "center"}
        ),
        html.P(id="molecule-name", style={"textAlign": "center"})
    ])
], className="mb-3")

app.layout = dbc.Container([
    dcc.Store(id="molecule-list"),
    # Store current config as default
    dcc.Store(id="current-config", data=load_hyperparameters()),
    dcc.Store(id="ga-running", data=False),
    # Interval to update live tracking every 1 second also used for checking ga status
    dcc.Interval(id="interval-component",interval=1000,n_intervals=0),
    header,
    dbc.Row([
        # Left columns for configuration panels
        dbc.Col([
            # Left config panel
            general_config_panel,
            dataset_toggle_panel,
            fitness_panel,
        ], width=2),
        dbc.Col([
            # Right config panel
            selection_panel,
            mutation_panel,
            crossover_panel,
            control_buttons_panel,
        ], width=2),
        # Right side
        dbc.Col([
            dbc.Row([
                dbc.Col(live_tracking_panel)
            ]),
            dbc.Row([
                dbc.Col(molecule_image_section, width=4),
                dbc.Col(results_panel, width=8),
            ])
        ], width=8)
    ])
], fluid=True)

@app.callback(
    Output("current-config", "data"),
    [
        Input("pop-size", "value"),
        Input("num-gen", "value"),
        Input("carrying-capacity", "value"),
        Input("mutation-weights-store", "data"),
        Input("selection-method", "value"),
        Input("crossover-methods", "value"),
        Input("elitism-weight", "value"),
        Input("data-toggle", "value"),
    ],
    State("current-config", "data")
)
def update_config(population_size, num_generations, carrying_capacity, mutation_weights,
                  selection_method, crossover_methods, elitism_weight, data_source, 
                  current_config):
    # Slider updates result in the config converting to a string
    # So, it must be manually converted back to dict here
    current_config = {}
    current_config["population_size"] = population_size
    current_config["num_generations"] = num_generations
    current_config["carrying_capacity"] = carrying_capacity
    current_config["mutation_weights"] = mutation_weights
    current_config["selection_method"] = selection_method
    current_config["crossover_methods"] = crossover_methods
    current_config["elitism_weight"] = elitism_weight
    current_config["data_source"] = data_source

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
    if not n_clicks or is_running:
        return is_running

    set_ga_active(True)
    # Open a new thread to runt the GA
    def run_ga_thread():
        config = load_config()
        ga = GeneticAlgorithm(config)
        ga.start_ga()
            
    thread = Thread(target=run_ga_thread)
    thread.start()
    return True

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
                x=0.5, y=0.5,
                xref="paper", yref="paper",
                font=dict(size=30, color="rgba(153, 179, 255,0.2)"),
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
        tickfont=dict(size=14, color="lightblue", family="Arial Bold"),
        range=[0,1],
        fixedrange=True,
        showgrid=False,
        zeroline=False,
    )
    return last_diversity_fig

@app.callback(
        Output("graph-3", "figure", allow_duplicate=True),
        Input("interval-component", "n_intervals"),
        State("ga-running", "data"),
        prevent_initial_call=True
)
def update_fitness_chart(n_iternvals, is_running):
    fitness_log = get_fitness_history_log()

    if not is_ga_active or not fitness_log:
        return px.line(title="Fitness Metrics")
    
    generations = list(range(1, len(fitness_log) + 1))
    min_fitness = [fit[0] for fit in fitness_log]
    mean_fitness = [fit[1] for fit in fitness_log]
    max_fitness = [fit[2] for fit in fitness_log]

    fig = px.line(labels={"x": "", "y": ""})
    fig.add_scatter(x=generations, y=max_fitness, mode="lines+markers", name="Max Fitness")
    fig.add_scatter(x=generations, y=mean_fitness, mode="lines+markers", name="Mean Fitness")
    fig.add_scatter(x=generations, y=min_fitness, mode="lines+markers", name="Min Fitness")
    # Reused format from Population Diversity chart.
    fig.update_layout(
        annotations=[
            dict(
                text="Fitness Progress",
                x=0.5, y=0.5,
                xref="paper", yref="paper",
                font=dict(size=30, color="rgba(153, 179, 255,0.2)"),
                showarrow=False
            ),
        ],
        legend=dict(
            orientation="h",
            y=0,
            x=0.5,
            xanchor="center",
            yanchor="top",
            font=dict(color="white"),
            bgcolor="rgba(0,0,0,0)"
        ),
        autosize=True,
        margin=dict(l=0,r=0,t=0,b=0,pad=0),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)"
    )
    fig.update_xaxes(
        showticklabels=False,
        ticks="",
        fixedrange=True,  # Prevents zooming
        showgrid=False,
        zeroline=False,
    )
    fig.update_yaxes(
        showticklabels=True,
        ticks="outside",
        tickvals=[0.0,0.5,1.0],
        tickfont=dict(size=14, color="lightblue", family="Arial Bold"),
        range=[0,1],
        fixedrange=True,
        showgrid=False,
        zeroline=False,
    )
    return fig

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
            df = evaluate_mols(population, top_n=50, reference_set=reference_set)
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

        mol = Chem.MolFromSmiles(smiles)

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

@app.callback(
    Output("mutation-percentage-display", "children"),
    Output("crossover-percentage-display", "children"),
    Output("avg-mutation-percentage-display", "children"),
    Output("avg-crossover-percentage-display", "children"),
    Input("interval-component", "n_intervals")
)
def update_success_rates(n_intervals):
    mutation_rate = get_latest_mutation_rates()
    crossover_rate = get_latest_crossover_rates()

    mutation_log = get_mutation_log()
    crossover_log = get_crossover_log()
    avg_mutation_rate = round(sum(mutation_log) / len(mutation_log) if mutation_log else 0)
    avg_crossover_rate = round(sum(crossover_log) / len(crossover_log) if crossover_log else 0)

    mutation_percentage_text = f"{mutation_rate:.0f}%"
    crossover_percentage_text = f"{crossover_rate:.0f}%"
    avg_mutation_percentage_text = f"{avg_mutation_rate:.0f}%"
    avg_crossover_percentage_text = f"{avg_crossover_rate:.0f}%"


    return mutation_percentage_text, crossover_percentage_text, avg_mutation_percentage_text, avg_crossover_percentage_text

@app.callback(
        Output("current-generation-number-display", "children"),
        Output("current-population-size-display", "children"),
        Input("interval-component", "n_intervals")
)
def update_ga_stats(n_intervals):
    current_generation_number = get_current_generation_number()
    current_population_size = get_current_population_size()

    current_generation_number_text = f"{current_generation_number:.0f}"
    current_population_size_text = f"{current_population_size:.0f}"
    
    return current_generation_number_text, current_population_size_text

@app.callback(
    Output("fitness-function-selection", "value"),
    Input("fitness-function-selection", "value")
)
def update_selected_fitness(selected_functions):
    update_selected_fitness_functions(selected_functions)
    return selected_functions

@app.callback(
    Output("mutation-weights-store", "data"),
    [
        Input("atomic-input", "value"),
        Input("functional-group-input", "value"),
        Input("ring-input", "value"),
        Input("fragment-input", "value")
    ]
)
def update_mutation_probabilities_callback(atomic, functional, ring, fragment):
    """
    Updates global mutation probabilites dictionary dynamically based on user input
    """
    updated_weights = {
        "atomic_substitution": atomic,
        "functional_group": functional,
        "ring": ring,
        "framgment": fragment,
    }

    update_mutation_probabilities(updated_weights)
    return updated_weights

@app.callback(
    Output("basic-fitness-collapse", "is_open"),
    Output("advanced-fitness-collapse", "is_open"),
    Input("fitness-advanced-toggle", "value")
)
def toggle_fitness_mode(toggle_val):
    if "advanced" in toggle_val:
        return False, True
    return True, False

@app.callback(
    Output("fitness-advanced-toggle", "style"),
    Input("fitness-advanced-toggle", "value"),
    Input("qed-tuning", "value"),
    Input("sa-tuning", "value"),
    Input("mw-tuning", "value"),
    Input("mw-filter-min", "value"),
    Input("mw-filter-max", "value"),
    Input("logp-filter-min", "value"),
    Input("logp-filter-max", "value"),
    Input("filter-selection", "value"),
)
def update_fitness_settings(advanced_toggle, qed, sa, mw, mw_min, mw_max, logp_min, logp_max, filters):
    # Update Weights
    tuning_weights = {
        "qed": qed,
        "sa": sa,
        "mol_weight": mw
    }
    update_tuning_weights(tuning_weights)

    if "advanced" in advanced_toggle:
        # Update Filters
        update_active_filters(filters)
        update_mw_range([mw_min, mw_max])
        update_logp_range([logp_min, logp_max])
    else:
        update_active_filters([])

    return dash.no_update

if __name__ == "__main__":
    app.run_server(debug=True)