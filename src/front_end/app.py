from flask import Flask, request, redirect, url_for, render_template
from controller import start_ga
from threading import Thread
import yaml
import time

app = Flask(__name__)

@app.route("/", methods=['GET', 'POST'])
def index():
    if request.methods == "POST":
        hyperparameters = request.form.to_dict()
        save_hyperparameters(hyperparameters)
        return redirect(url_for("index"))
    
    hyperparameters = load_hyperparameters()
    return render_template(index.html, hyperparameters=hyperparameters)

def load_hyperparameters():
    with open("../../config.yaml","r") as file:
        return yaml.safe_load(file)

def save_hyperparameters(hyperparameters):
    with open("../../config.yaml", "w") as file:
        yaml.dump(hyperparameters, file)
        
is_running = False

@app.route("/start")
def start_ga():
    global is_running
    if not is_running:
        is_running = True
        thread = Thread(target=run_ga)
        thread.start()
        return "GA Started"
    else:
        return "GA already running"

def stop_ga():
    global is_running
    is_running = False
    return "GA manually stopped"
        
        
def run_ga():
    global is_running
    is_running = True

    start_ga()
        
if __name__ == "__main__":
    app.run(debug=True)