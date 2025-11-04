import pandas as pd

df=pd.read_csv("final_merged_dataset_fuzzy_ATC_fullDB.csv",sep=",")

# Remplacer les NaN dans 'ATC' par chaîne vide
df['ATC'] = df['ATC'].fillna('')

# Créer une colonne 'ATC_list' avec une liste des codes séparés par ';'
df['ATC_list'] = df['ATC'].apply(lambda x: [code for code in x.split(';') if code])

# Remplacer la colonne 'ATC' par une chaîne formée à partir de cette liste (codes séparés par ';')
df['ATC'] = df['ATC_list'].apply(lambda x: ';'.join(x))

# Supprimer la colonne temporaire 'ATC_list'
df = df.drop(columns=['ATC_list'])


# Afficher le DataFrame pour vérification
print(df)
df.to_csv("df.csv", index=False)
print(df.info())