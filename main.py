from dash import Dash, html, dcc, Input, Output, State
import pandas as pd

user_list = list(range(183)) #list(set(import_df["user"])) # List of all users 
transportation_list = ['train', 'taxi', 'walk', 'bus', 'subway', 'airplane', 'car', 'bike', 'motorcycle', 'run', 'boat'] # List of all transportation's
graphing_method_list = ["Angle"] # List of all graphing methods TBI
presets_list = ["very_short", "short", "medium", "car", "walk", "PT"] # List of all presets

# Creates the dash app
app = Dash()

app.layout = html.Div([
    dcc.RadioItems(options=["Presets", "Custom"], value="Presets", id="data_type_radio", style={"display": "inline-block"}, inline=True),

    html.Div([html.Div("Preset:", style={"display": "inline-block"}),
    dcc.Dropdown(options=presets_list, value=presets_list[0], id="preset_dropdown", 
                 style={"display": "inline-block"})],
                 id = "presets_div", style={"display": "inline-block", "width": "20%"}),

    # Filter on user if filtering
    html.Div([html.Div("User:", style={"display": "inline-block"}),
    dcc.Dropdown(options=["All"] + user_list, value="All", id="user_filter_dropdown", 
                 style={"display": "inline-block", "width": "45%"}),
    
    # Filter on transportation if filtering
    html.Div("Transportation:", style={"display": "inline-block"}),
    dcc.Dropdown(options=["All"] + transportation_list, value="All", id="transportation_filter_dropdown", 
                 style={"display": "inline-block", "width": "45%"}),], 
                 id = "filtering_div", style = {"display": "none", "width": "40%"}),
    
    # Choose the method for graphing
    html.Div("Graphing Method:", style={"display": "inline-block"}),
    dcc.Dropdown(options=graphing_method_list, value=graphing_method_list[0], id="graphing_method_dropdown", 
                 style={"display": "inline-block", "width": "20%"}),
    
    # Runs the application
    html.Button("Run", id="run_button", style={"display": "inline-block"}),

    # Space for the two graphs
    dcc.Graph(id="base_visualization_graph", config={"displayModeBar": False}, style={"display": "inline-block", "width": "40%"}),
    dcc.Graph(id="graph_visualization_graph", config={"displayModeBar": False}, style={"display": "inline-block", "width": "40%"}),
])

@app.callback(
    Output("filtering_div", "style"),
    Output("presets_div", "style"),
    Input("data_type_radio", "value"),
)
def active_filter(data_type):
    if data_type == "Presets":
        return {"display": "none"}, {"display": "inline-block"}
    elif data_type == "Custom":
        return {"display": "inline-block"}, {"display": "none"}
    
@app.callback(
    Output("base_visualization_graph", "figure"),
    Output("graph_visualization_graph", "figure"),
    Input("run_button", "n_clicks"),
    State("data_type_radio", "value"),
    State("preset_dropdown", "value"),
    State("user_filter_dropdown", "value"),
    State("transportation_filter_dropdown", "value"),
    prevent_initial_call=True
)
def run(dummy, data_type, preset, user, transportation):
    if data_type == "Presets":
        df = pd.read_csv(f"presets/{preset}")

    elif data_type == "Custom":
        df = pd.read_csv("presets/full_data.csv")

        if user != "All":
            df = df[df["user"] == user]
        
        if transportation != "All":
            df = df[df["transportation"] == transportation]

    # Create the base visualization using the filtered 'df' dataframe TBI
    base_visualization = None 

    # Create the graphed visualization using the filtered 'df' dataframe and the graphing method TBI
    graphing_visualization = None 

    return base_visualization, graphing_visualization


if __name__ == "__main__":
    # Runs the app
    app.run_server(debug=True)