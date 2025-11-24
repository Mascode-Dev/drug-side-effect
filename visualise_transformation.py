import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition import PCA

# 1. Chargement des données
# Remplacez par le chemin de votre fichier si nécessaire
df = pd.read_csv('cleaned_drug_side_effects.csv')

# ---------------------------------------------------------
# PARTIE 1 : La "Carte" des médicaments (PCA)
# ---------------------------------------------------------

# A. Encodage TF-IDF
# On demande des n-grams de 3 à 5 lettres
tfidf = TfidfVectorizer(analyzer='char', ngram_range=(3, 5))
X_names = tfidf.fit_transform(df['drug_name'])

# B. Réduction de dimension (PCA)
# On passe de milliers de dimensions à 2 (x et y) pour pouvoir dessiner
pca = PCA(n_components=2)
coords = pca.fit_transform(X_names.toarray())

# Création d'un DataFrame temporaire pour le graphique
plot_df = pd.DataFrame(coords, columns=['PC1', 'PC2'])
plot_df['drug_name'] = df['drug_name']

# C. Création du graphique
plt.figure(figsize=(14, 9))
# On dessine d'abord tous les points en gris pour le fond
sns.scatterplot(data=plot_df, x='PC1', y='PC2', alpha=0.5, s=40, color='lightgrey')

# On surligne quelques familles spécifiques pour montrer l'efficacité
familles_a_voir = ['azole', 'cillin', 'pril'] # Ex: antifongiques, pénicillines, hypertenseurs
couleurs = ['red', 'green', 'orange']

for suffixe, couleur in zip(familles_a_voir, couleurs):
    # On filtre les médicaments qui contiennent ce suffixe
    mask = plot_df['drug_name'].str.contains(suffixe)
    subset = plot_df[mask]
    
    if not subset.empty:
        # On les dessine en couleur par-dessus
        plt.scatter(subset['PC1'], subset['PC2'], color=couleur, label=f"Contient '{suffixe}'", s=60, edgecolor='black')
        
        # Optionnel : Ajouter le texte pour quelques points seulement (pour ne pas surcharger)
        for i in range(min(5, len(subset))): # Max 5 labels par famille
            row = subset.iloc[i]
            plt.text(row['PC1'], row['PC2']+0.005, row['drug_name'], fontsize=9, color=couleur, weight='bold')

plt.title("Visualisation de l'encodage TF-IDF : Regroupement automatique par similarité de nom", fontsize=14)
plt.xlabel("Composante Principale 1")
plt.ylabel("Composante Principale 2")
plt.legend()
plt.grid(True, linestyle='--', alpha=0.3)
plt.show()


# ---------------------------------------------------------
# PARTIE 2 : Le Zoom "Heatmap"
# ---------------------------------------------------------

# A. Sélection d'un petit échantillon pour la démo
# On prend 5 médicaments 'azole' et on ajoute 'aspirin' pour comparer
target_drugs = df[df['drug_name'].str.contains('azole')].head(5)['drug_name'].tolist()
target_drugs.append('aspirin') 

# B. On refait un TF-IDF local juste sur ces mots pour avoir une matrice lisible
tfidf_subset = TfidfVectorizer(analyzer='char', ngram_range=(3, 5))
X_subset = tfidf_subset.fit_transform(target_drugs)
feature_names = tfidf_subset.get_feature_names_out()

# C. Conversion en DataFrame
df_heatmap = pd.DataFrame(X_subset.toarray(), index=target_drugs, columns=feature_names)

# D. Filtrage pour l'affichage
# On ne garde que les 20 n-grams les plus forts (sinon le graphique est trop large)
top_cols = df_heatmap.max().sort_values(ascending=False).head(20).index
df_heatmap_limited = df_heatmap[top_cols]

# E. Affichage Heatmap
plt.figure(figsize=(12, 6))
sns.heatmap(df_heatmap_limited, cmap='viridis', annot=True, fmt=".2f", cbar_kws={'label': 'Score TF-IDF'})
plt.title("Zoom : Pourquoi l'algorithme regroupe ces médicaments ? (Analyse des n-grams communs)")
plt.xlabel("Segments de caractères (features)")
plt.ylabel("Médicament")
plt.tight_layout()
plt.show()