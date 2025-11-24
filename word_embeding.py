import pandas as pd
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.feature_extraction.text import TfidfVectorizer

# 1. Chargement
df = pd.read_csv('cleaned_drug_side_effects.csv')

# Fonction utilitaire pour transformer les chaînes "A;B;C" en listes ['A', 'B', 'C']
def split_items(x):
    return x.split(';') if pd.notnull(x) else []

# --- A. Encodage des Targets (Variables explicatives) ---
# On transforme la colonne en liste
df['targets_list'] = df['targets'].apply(split_items)

# MultiLabelBinarizer crée une colonne par cible possible
mlb_targets = MultiLabelBinarizer()
X_targets = mlb_targets.fit_transform(df['targets_list'])
print(f"Shape des features 'targets': {X_targets.shape}") 
# Résultat attendu : (1052, ~1148)

# --- B. Encodage des Noms (Variables explicatives) ---
# On utilise TF-IDF sur les caractères (n-grams de taille 3 à 5)
# Capture des racines comme "meth", "oxy", "zole"
# analyzer='char' : Au lieu de lire des mots entiers (comme "chat", "chien"), l'algorithme va lire des caractères (lettres). C'est crucial car les noms de médicaments n'ont pas de "phrases", 
# le sens est caché dans le mot lui-même. ngram_range=(3, 5) : Cela demande à l'algorithme de créer une fenêtre glissante qui prend des groupes de 3, 4 et 5 lettres consécutives.
tfidf = TfidfVectorizer(analyzer='char', ngram_range=(3, 5))

#Fit (Apprendre) : Il parcourt toute votre colonne drug_name, découpe tous les noms en n-grams (3 à 5 lettres) et construit un gigantesque dictionnaire de toutes les combinaisons existantes (ex: azep, zepa, epam...)
#Transform (Encoder) : Il remplace chaque nom de médicament par un vecteur numérique (une suite de chiffres).
X_names = tfidf.fit_transform(df['drug_name'])

print(f"Shape des features 'names': {X_names.shape}")

# --- C. Encodage des Effets Secondaires (Cible à prédire) ---
df['side_effects_list'] = df['side_effect'].apply(split_items)

mlb_effects = MultiLabelBinarizer()
y = mlb_effects.fit_transform(df['side_effects_list'])
print(f"Shape de la cible 'y': {y.shape}")
# Résultat attendu : (1052, ~5735)

# --- D. Fusion pour le modèle ---
# Pour obtenir votre matrice X finale, vous concaténez les features
import numpy as np
X_final = np.hstack([X_targets, X_names.toarray()])

print(f"Matrice finale pour entraînement : {X_final.shape}")
print(f" Premieres lignes de la matrice finale X_final :\n{X_final[:5]}")  # Affiche les 5 premières lignes et 5 premières colonnes