import streamlit as st
import array
import math
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem
from chembl_webresource_client.new_client import new_client
from tqdm.auto import tqdm
import os

HERE = Path(os.getcwd())
DATA = HERE / "data"
 
st.title("ChemBL Bioactivity Data App")

# User Input for Target
uniprot_id = st.text_input("Enter Uniprot ID for the target (e.g., P00533):", "P00533")

# Get Target Data
targets_api = new_client.target
targets = targets_api.get(target_components__accession=uniprot_id).only(
    "target_chembl_id", "organism", "pref_name", "target_type"
)
targets = pd.DataFrame.from_records(targets)
if not targets.empty:
    st.subheader("Target Information:")
    st.write(targets)
 
    target = targets.iloc[0]
    chembl_id = target.target_chembl_id
    st.progress(15)
    st.subheader("fetching bioactivity data for target...")

    # Get Bioactivity Data
    bioactivities_api = new_client.activity
    bioactivities = bioactivities_api.filter(
        target_chembl_id=chembl_id, tye="IC50", relation="=", assay_type="B"
    ).only(
        "activity_id",
        "assay_chembl_id",
        "assay_description",
        "assay_type",
        "molecule_chembl_id",
        "type",
        "standard_units",
        "relation",
        "standard_value",
        "target_chemnl_id",
        "target_organism",
    )

    bioactivities_df = pd.DataFrame.from_dict(bioactivities)
    # Preprocess Bioactivity Data
    st.progress(30)
    st.subheader("Raw Bioactivity Data:")
    st.write(bioactivities_df.head())
    bioactivities_df.drop(['units','value'], axis=1, inplace=True)
    bioactivities_df = bioactivities_df.astype({"standard_value": 'float64'})
    bioactivities_df.dropna(axis=0, how='any', inplace=True)
    bioactivities_df = bioactivities_df[bioactivities_df["standard_units"] == "nM"]
    bioactivities_df.drop_duplicates("molecule_chembl_id", keep="first", inplace=True)
    bioactivities_df.reset_index(drop=True, inplace=True)
    bioactivities_df.rename(
    columns={
        "standard_value": "IC50",
        "standard_units": "units",
    },
    inplace=True,
)

    st.progress(40)
    st.subheader("Processed Bioactivity Data:")
    st.write(bioactivities_df.head())

    # Fetch Compound Data
    compounds_api = new_client.molecule
    compounds_provider = compounds_api.filter(
        molecule_chembl_id__in=list(bioactivities_df["molecule_chembl_id"])
    ).only("molecule_chembl_id", "molecule_structures")
    compounds = list(tqdm(compounds_provider))
    compounds_df = pd.DataFrame.from_records(compounds)

    # Preprocess Compound Data
    st.progress(65)
    #Remove entries with missing entries
    compounds_df.dropna(axis=0, how='any', inplace=True)

    #Delete duplicate molecules (from mol chembl ID)
    compounds_df.drop_duplicates("molecule_chembl_id", keep="first", inplace=True)

    #now we only select smiles
    canonical_smiles = []
    for i, compounds in compounds_df.iterrows():
        try:
            canonical_smiles.append(compounds["molecule_structures"]["canonical_smiles"])
        except KeyError:
            canonical_smiles.append(None)
    compounds_df["smiles"] = canonical_smiles
    compounds_df.drop("molecule_structures", axis=1, inplace=True)

    #print(f'DataFrame shape: {compounds_df.shape}')
    #just to be sure we will also remove all data without a smiles string
    compounds_df.dropna(axis=0, how='any', inplace=True)




    st.subheader("Processed Compound Data:")
    st.write(compounds_df.head())
 
    # Merge DataFrames
    st.progress(90)
    output_df = pd.merge(
        bioactivities_df[["molecule_chembl_id", "IC50", "units"]],
        compounds_df,
        on="molecule_chembl_id",
    )

    # Convert IC50 to pIC50
    output_df["pIC50"] = output_df["IC50"].apply(lambda x: 9 - math.log10(x))

    st.subheader("Final Merged Dataset:")
    st.write(output_df.head())

    # Draw Histogram
    st.progress(100)
    st.subheader("Histogram of pIC50 Values:")
    plt.hist(output_df["pIC50"])
    plt.xlabel("pIC50")
    plt.ylabel("Frequency")
    plt.title("Histogram of pIC50 Values")
    st.pyplot(plt)

 
    # Sort by pIC50 and save to CSV
    st.subheader("Data preview:")
    output_df_sorted = output_df.sort_values(by="pIC50", ascending=False)

    # Add molecule column
    PandasTools.AddMoleculeColumnToFrame(output_df_sorted, smilesCol="smiles")

    # Reset index
    output_df_sorted.reset_index(drop=True, inplace=True)

   

    # Prepare saving the dataset: Drop the ROMol column
    output_df_sorted = output_df_sorted.drop("ROMol", axis=1)
    st.write(output_df_sorted.head())
    st.subheader("Bar Chart of pIC50 Values for the First 5 Molecules:")
    chart_data = output_df.head(5)
    st.bar_chart(chart_data.set_index('molecule_chembl_id')['pIC50'])
    st.subheader("Download Data:")
    st.dataframe(output_df_sorted)
    st.success("Dataset processed successfully!")
else:
    st.warning("No target information available for the provided Uniprot ID.")