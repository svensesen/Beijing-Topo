from dash import Dash, html, dcc, Input, Output
import pandas as pd

# Counter to be able to send signals
dataset_filter_counter = 0

import_df = pd.DataFrame() # Importing the dataframe TBI
user_list = ["000", "001", "002"] # List of all users TBI (should get this from the df)
transportation_list = ["Walk", "Bike", "Bus"] # List of all transportation"s TBI (should get this from the df)
graphing_method_list = ["Angle"] # List of all graphing methods TBI

# Creates the dash app
app = Dash()

app.layout = html.Div([
    html.Div("User:", style={"display": "inline-block"}),
    dcc.Dropdown(options=["All"] + user_list, value="All", id="user_filter_dropdown", style={"display": "inline-block", "width": "20%"}),
    html.Div("Transportation:", style={"display": "inline-block"}),
    dcc.Dropdown(options=["All"] + transportation_list, value="All", id="transportation_filter_dropdown", style={"display": "inline-block", "width": "20%"}),
    html.Div("Graphing Method:", style={"display": "inline-block"}),
    dcc.Dropdown(options=graphing_method_list, value=graphing_method_list[0], id="graphing_method_dropdown", style={"display": "inline-block", "width": "20%"}),

    dcc.Graph(id="base_visualization_graph", config={"displayModeBar": False}, style={"display": "inline-block", "width": "40%"}),
    dcc.Graph(id="graph_visualization_graph", config={"displayModeBar": False}, style={"display": "inline-block", "width": "40%"}),

    dcc.Store(data=[0], id="dataset_filter"), # Keeps track of whenever the dataset is filtered
])

@app.callback(
    Output("dataset_filter", "data"),
    Input("user_filter_dropdown", "value"),
    Input("transportation_filter_dropdown", "value"),
)
def filter(user, transportation):
    global dataset_filter_counter, df
    dataset_filter_counter += 1

    # Use user and transportation to filter TBI
    df = import_df

    return [dataset_filter_counter]


@app.callback(
    Output("base_visualization_graph", "figure"),
    Output("graph_visualization_graph", "figure"),
    Input("graphing_method_dropdown", "value"),
    Input("dataset_filter", "data"),
)
def filter(graphing_method, dummy):
    # Create the base visualization using the filtered 'df' dataframe TBI
    base_visualization = None 

    # Create the graphed visualization using the filtered 'df' dataframe and the graphing method TBI
    graphing_visualization = None 

    return base_visualization, graphing_visualization


if __name__ == "__main__":
    # Runs the app
    app.run_server(debug=True)