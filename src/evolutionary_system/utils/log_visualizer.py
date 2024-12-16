import os
import pandas as pd
import matplotlib.pyplot as plt

def read_log():
    cd = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(cd, "metrics_log.csv")
    
    try:
        df = pd.read_csv(file_path)
        return df
    except FileNotFoundError:
        print(f"Error - {file_path} not found.")
        return None
    
def plot_runtime(df):
    try:
        plt.figure(figsize=(10,6))
        plt.plot(df["run_id"], df["runtime"], color='g', linewidth=2)
        plt.xlabel("Run ID")
        plt.ylabel("Runtime (s)")
        plt.title("Runtime")
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.show()
    except Exception as e:
        print(f"Error plotting runtime: {e}")

def main():
    df = read_log()
    if df is not None:
        plot_runtime(df)

if __name__ == "__main__":
    main()
