# 🔮 Drug Side Effect Predictor

This project aims to predict the potential side effects of a drug molecule based on its **name (encoded chemical structure)** and its **biological targets**. It utilizes Natural Language Processing (NLP) techniques and multi-label supervised learning.

## 📋 Table of Contents

1. Project Overview
2. Repository Structure
3. Installation
4. Data Pipeline
5. Modeling and Performance
6. Dashboard Usage
7. Contribution of each member

## 🚀 Project Overview

The project merges data from **STITCH**, **SIDER (MedDRA)**, and **DrugBank** to create a comprehensive dataset. This dataset associates:

* Drug names.
* Anatomical Therapeutic Chemical (ATC) codes.
* Biological targets.
* Indications and side effects.

An optimized **Logistic Regression (OneVsRest)** model is trained to predict the presence or absence of the 300 most frequent side effects.

## 📂 Repository Structure

* `Fusionetconcat.py`: Script for initial merging of STITCH and DrugBank databases, featuring XML parsing and Fuzzy Matching.
* `Firstcleaning.py`: Script for cleaning and formatting ATC codes.
* `data_cleaning.ipynb`: Exploratory data analysis and final dataset cleaning.
* `word_embedding.ipynb`: Data encoding (TF-IDF on character n-grams) and model training/comparison.
* `dashboard.py`: Interactive Dash application to test the model in real-time.
* `models_export/`: Directory containing trained models and vectorizers in `.pkl` format.

## 🛠 Installation

1. Clone the repository.
2. Install the necessary dependencies:

```bash
pip install pandas numpy scikit-learn dash joblib matplotlib plotly rapidfuzz

```

3. Ensure you have the source data files or the generated `cleaned_drug_side_effects.csv`.

## 📊 Data Pipeline

### 1. Preprocessing

* **XML Parsing**: Extraction of targets and descriptions from the DrugBank database.
* **Fuzzy Matching**: Reconciling drug names between STITCH and DrugBank using the `token_sort_ratio` algorithm.
* **Cleaning**: Removal of missing values and filtering of non-exploitable columns such as IDs and long descriptions.

### 2. Encoding (Embedding)

* **Biological Targets**: Binary encoding via `MultiLabelBinarizer`.
* **Molecule Names**: Use of **TF-IDF on characters** (n-grams of 3 to 5 letters) to capture chemical roots like *-zole* or *-meth*.

## 🤖 Modeling and Performance

The project compared several approaches:

* **Baseline (Logistic Regression)**: Optimized via `GridSearchCV`.
* **Random Forest**.
* **Voting Classifier**.

**Result:** The Logistic Regression model achieved the best performance with a **Micro F1-Score of approximately 0.49** on the top 300 side effects. This is highly significant for a multi-label problem of this complexity.

## 🖥 Dashboard Usage

To launch the prediction interface:

```bash
python dashboard.py

```

Once launched, access `http://127.0.0.1:8050/` in your web browser. You can:

1. Enter a molecule name (e.g., Aspirin).
2. Select one or more biological targets from the dropdown menu.
3. Instantly receive a list of potential side effects detected by the model.

---

## 👥 Contibutions of each member
Ghita Fadili : Dataset selection and merging, Scientific Papers & References\
Thomas Jego : Model training, Ensemble models, Voting, Model Comparison\
Victoire Dezombre : Data Exploration & Data Cleaning\
Léandre durand-Terasson : Word Embedding (TF-IDF, syntactic similarity, MultiLabelBinarizer)
