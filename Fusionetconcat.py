import pandas as pd
import xml.etree.ElementTree as ET
from rapidfuzz import process, fuzz

# -----------------------------
# 1️⃣ Charger STITCH / MEDDRA
# -----------------------------
drug_names = pd.read_csv("drug_names.tsv", sep="\t", header=None, names=["stitch_id", "drug_name"])
drug_atc = pd.read_csv("drug_ATC.tsv", sep="\t", header=None, names=["stitch_id", "ATC_code"])
drug_atc.drop_duplicates(subset=["stitch_id"], inplace=True)

meddra_se = pd.read_csv(
    "meddra_all_se.tsv", sep="\t", header=None,
    names=["stitch_id", "cid2", "umls_id", "concept_type", "meddra_id", "side_effect"]
)
meddra_ind = pd.read_csv(
    "meddra_all_indications.tsv", sep="\t", header=None,
    names=["stitch_id", "umls_id", "source", "disease_name", "concept_type", "meddra_id", "disease_term"]
)

# -----------------------------
# 2️⃣ Agréger indications et effets secondaires
# -----------------------------
ind_agg = meddra_ind.groupby("stitch_id")["disease_term"].apply(lambda x: ";".join(sorted(set(x.dropna())))).reset_index()
se_agg = meddra_se.groupby("stitch_id")["side_effect"].apply(lambda x: ";".join(sorted(set(x.dropna())))).reset_index()

# -----------------------------
# 3️⃣ Fusion STITCH
# -----------------------------
merged_stitch = (
    drug_names.merge(drug_atc, on="stitch_id", how="left")
    .merge(ind_agg, on="stitch_id", how="left")
    .merge(se_agg, on="stitch_id", how="left")
)
merged_stitch['drug_name_lower'] = merged_stitch['drug_name'].str.lower().str.strip()

# -----------------------------
# 4️⃣ Parser DrugBank XML COMPLET
# -----------------------------
xml_path = "full database.xml"
ns = {"db": "http://www.drugbank.ca"}

context = ET.iterparse(xml_path, events=("start", "end"))
_, root = next(context)

data = []
count = 0
for event, elem in context:
    if event == "end" and elem.tag == "{http://www.drugbank.ca}drug":
        drug_type = elem.attrib.get("type", "")
        if drug_type not in ["small molecule", "biotech"]:
            root.clear()
            continue

        ids = [d.text for d in elem.findall("db:drugbank-id", namespaces=ns) if d.text]
        primary_id = elem.findtext("db:drugbank-id[@primary='true']", namespaces=ns)
        drugbank_id = primary_id if primary_id else (ids[0] if ids else None)

        name = elem.findtext("db:name", namespaces=ns)
        description = elem.findtext("db:description", namespaces=ns)
        groups = [g.text for g in elem.findall(".//db:groups/db:group", namespaces=ns) if g.text]
        atc_codes = [a.attrib.get("code") for a in elem.findall(".//db:atc-code", namespaces=ns)]
        synonyms = [s.text for s in elem.findall(".//db:synonyms/db:synonym", namespaces=ns) if s.text]
        targets = [t.findtext("db:name", namespaces=ns) for t in elem.findall(".//db:targets/db:target", namespaces=ns)]

        data.append({
            "drugbank_id": drugbank_id,
            "name": name,
            "description": description,
            "type": drug_type,
            "groups": "; ".join(groups),
            "atc_codes": "; ".join(atc_codes),
            "synonyms": "; ".join(synonyms),
            "targets": "; ".join(targets)
        })

        count += 1
        if count % 1000 == 0:
            print(f"Parsed {count} drugs...")

        # vider le root pour la mémoire
        root.clear()

df_drugbank = pd.DataFrame(data)
df_drugbank['name_lower'] = df_drugbank['name'].str.lower().str.strip()

print(f"✅ Total DrugBank drugs parsed: {len(df_drugbank)}")

# -----------------------------
# 5️⃣ Fuzzy Matching STITCH ↔ DrugBank
# -----------------------------
drugbank_names = df_drugbank['name_lower'].tolist()

def get_best_match(name, choices, threshold=80):
    match = process.extractOne(name, choices, scorer=fuzz.token_sort_ratio)
    if match and match[1] >= threshold:
        return match[0]
    return None

merged_stitch['matched_name'] = merged_stitch['drug_name_lower'].apply(lambda x: get_best_match(x, drugbank_names))

# -----------------------------
# 6️⃣ Fusion finale
# -----------------------------
final_df = merged_stitch.merge(
    df_drugbank,
    left_on='matched_name',
    right_on='name_lower',
    how='left'
)

# -----------------------------
# 7️⃣ Fusion ATC en une seule colonne
# -----------------------------
def merge_atc(row):
    codes = []
    if pd.notna(row['ATC_code']):
        codes.append(row['ATC_code'])
    if pd.notna(row['atc_codes']):
        codes.append(row['atc_codes'])
    if codes:
        codes_unique = sorted(set(";".join(codes).split(";")))
        return ";".join(codes_unique)
    return None

final_df['ATC'] = final_df.apply(merge_atc, axis=1)

# -----------------------------
# 8️⃣ Nettoyage
# -----------------------------
cols_to_drop = ['ATC_code', 'atc_codes', 'drug_name_lower', 'name_lower', 'matched_name']
final_df.drop(columns=[c for c in cols_to_drop if c in final_df.columns], inplace=True)

# -----------------------------
# 9️⃣ Sauvegarde
# -----------------------------
final_df.to_csv("final_merged_dataset_fuzzy_ATC_fullDB.csv", index=False, encoding="utf-8")

print("✅ Fusion complète avec DrugBank complet et ATC combiné terminée !")
print(f"Dimensions finales : {final_df.shape}")
print(final_df.sample(5))
print(final_df.info())
