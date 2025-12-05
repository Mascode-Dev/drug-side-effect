import dash
from dash import dcc, html, Input, Output, State
import joblib
import numpy as np
import pandas as pd
import os

# =========================================================
# 1. CHARGEMENT DES MODÈLES ET OUTILS
# =========================================================
model_path = 'models_export'

try:
    print("Models loading...")
    model = joblib.load(os.path.join(model_path, 'best_model.pkl'))
    tfidf_name = joblib.load(os.path.join(model_path, 'tfidf_name.pkl'))
    mlb_targets = joblib.load(os.path.join(model_path, 'mlb_targets.pkl'))
    effects_names = joblib.load(os.path.join(model_path, 'filtered_effects_names.pkl'))
    
    print(f"Models loaded. The system knows {len(effects_names)} side effects.")
    
    available_targets = sorted(mlb_targets.classes_)

except FileNotFoundError as e:
    print(f"ERROR: {e}")
    print("Make sure you have exported 'filtered_effects_names.pkl' from the notebook.")
    available_targets = []
    model = None
    effects_names = []

# =========================================================
# 2. CONFIGURATION DASH
# =========================================================
app = dash.Dash(__name__, title="Drug Side Effects Predictor")

app.layout = html.Div(style={'fontFamily': 'Arial, sans-serif', 'maxWidth': '1000px', 'margin': '0 auto', 'padding': '20px'}, children=[
    html.Div(style={'textAlign': 'center', 'marginBottom': '40px'}, children=[
        html.H1("🔮 Prediction of Side Effects", style={'color': '#2c3e50'}),
        html.P("Enter the chemical composition (Name) and biological targets.", style={'color': '#7f8c8d'})
    ]),

    html.Div(style={'backgroundColor': '#f9f9f9', 'padding': '30px', 'borderRadius': '10px', 'boxShadow': '0 4px 6px rgba(0,0,0,0.1)'}, children=[
        html.Label("1. Name of the molecule / drug:", style={'fontWeight': 'bold'}),
        dcc.Input(id='input-name', type='text', placeholder='Ex: Aspirin...', style={'width': '100%', 'padding': '10px', 'marginBottom': '20px'}),

        html.Label("2. Biological Targets:", style={'fontWeight': 'bold'}),
        dcc.Dropdown(id='input-targets', options=[{'label': t, 'value': t} for t in available_targets], multi=True, placeholder="Select targets...", style={'marginBottom': '20px'}),
        html.Button('Run Prediction', id='predict-btn', n_clicks=0, style={'backgroundColor': '#3498db', 'color': 'white', 'padding': '15px', 'width': '100%', 'border': 'none', 'borderRadius': '5px', 'cursor': 'pointer'})
    ]),

    html.Div(id='result-container', style={'marginTop': '30px', 'display': 'none'}, children=[
        html.H3("Analysis Results:", style={'borderBottom': '2px solid #3498db'}),
        html.Div(id='output-prediction', style={'padding': '20px', 'backgroundColor': '#fff', 'border': '1px solid #eee'})
    ])
])

# =========================================================
# LOGIC FOR PREDICTION
# =========================================================
@app.callback(
    [Output('output-prediction', 'children'),
        Output('result-container', 'style')],
    [Input('predict-btn', 'n_clicks')],
    [State('input-name', 'value'),
        State('input-targets', 'value')]
)
def predict_side_effects(n_clicks, name_value, targets_value):
    if n_clicks == 0 or model is None:
        return "", {'display': 'none'}
    if not name_value:
        return html.Div("⚠️ Name missing.", style={'color': 'red'}), {'display': 'block'}

    try:
        # A. Preparation
        name_vector = tfidf_name.transform([name_value]).toarray()
        target_list = [targets_value] if targets_value else [[]]
        target_vector = mlb_targets.transform(target_list)
        features = np.hstack([target_vector, name_vector])

        # B. Prediction
        prediction_binary = model.predict(features) # Returns e.g.: [[0, 1, 0, 0, 1...]]

        # We will manually look where the "1"s are and take the corresponding name from our list.        
        indices_actifs = np.where(prediction_binary[0] == 1)[0]
        effects = [effects_names[i] for i in indices_actifs]
        
        # D. Display
        if len(effects) == 0:
            return html.Div([
                html.H4("✅ No side effects detected among the Top 500", style={'color': 'green'}),
                html.P("Note: The model only monitors the 500 most frequent effects.")
            ]), {'display': 'block'}
        
        return html.Div([
            html.H4(f"⚠️ {len(effects)} Side effects detected:", style={'color': '#e74c3c'}),
            html.Ul([html.Li(e) for e in effects])
        ]), {'display': 'block'}

    except Exception as e:
        return html.Div(f"Technical error: {str(e)}", style={'color': 'red'}), {'display': 'block'}

if __name__ == '__main__':
    app.run(debug=True, port=8050)