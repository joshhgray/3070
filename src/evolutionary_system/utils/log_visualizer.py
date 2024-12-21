import plotly.express as px
import pandas as pd

def plot_diversity(num_generations, diversity_log):
    '''
    Plot the diversity (Jaccard Index) for each generation
    '''
    df = pd.DataFrame({"Generation": num_generations, "Diversity": diversity_log})
    fig = px.line(df, x="Generation", y="Diversity", title="Diversity by Generation")
    fig.show()
