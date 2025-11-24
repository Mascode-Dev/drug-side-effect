import pandas as pd

# Load the dataset
file_path = 'cleaned_drug_side_effects.csv'
try:
    df = pd.read_csv(file_path)
    # Display the first few rows and info to understand the structure
    #print("First 5 rows:")
    #print(df.head().to_markdown(index=False, numalign="left", stralign="left"))
    print("\nDataFrame Info:")
    print(df.info())
    
    # Let's inspect the specific columns mentioned: drug_name, targets (composition?), side_effect
    #print("\nSample of 'targets' column:")
    #print(df['targets'].head(10).tolist())
    
    #print("\nSample of 'side_effect' column:")
    #print(df['side_effect'].head(10).tolist())
    
    print("\nSample of 'groups' column:")
    print(df['groups'].head(10).tolist())

except Exception as e:
    print(f"Error loading file: {e}")

# Function to count unique items in semicolon-separated columns
def count_unique_items(series):
    unique_items = set()
    for row in series.dropna():
        items = [item.strip() for item in row.split(';')]
        unique_items.update(items)
    return len(unique_items), list(unique_items)[:5]

n_targets, sample_targets = count_unique_items(df['targets'])
n_side_effects, sample_side_effects = count_unique_items(df['side_effect'])
n_drugs = df['drug_name'].nunique()

print("--------------------------------------------------------------------------------------------------------------")

print(f"Unique Targets: {n_targets}")
print(f"Sample Targets: {sample_targets}")
print(f"Unique Side Effects: {n_side_effects}")
print(f"Unique Drugs: {n_drugs}")

print("---------------------------------------------------------------------------------------------------------------")
# Check distribution of number of targets per drug
df['num_targets'] = df['targets'].apply(lambda x: len(str(x).split(';')) if pd.notnull(x) else 0)
print(df['num_targets'].describe())