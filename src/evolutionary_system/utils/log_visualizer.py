import plotly.express as px #type:ignore
import pandas as pd#type:ignore

def plot_diversity(num_generations, diversity_log):
    '''
    Plot the diversity (Jaccard Index) for each generation
    '''
    df = pd.DataFrame({"Generation": num_generations, "Diversity": diversity_log})
    fig = px.line(df, x="Generation", y="Diversity", title="Diversity by Generation")
    fig.show()
