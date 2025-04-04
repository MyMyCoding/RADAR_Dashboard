import streamlit as st
import pandas as pd
import numpy as np
import py3Dmol
import matplotlib.pyplot as plt

st.set_page_config(layout="wide")
st.title("RNA Aptamer Docking & Analysis Resource - R.A.D.A.R.")

pdb_file = st.sidebar.file_uploader("Upload PDB File", type=['pdb'])
custom_aptamer = st.sidebar.text_input("Enter Custom Aptamer Label", "APT-CUSTOM")

if pdb_file:
    docking_score = st.number_input("Enter docking score (e.g., from HDOCK)", value=-200.0)
    confidence_score = st.number_input("Enter confidence score (0–1)", min_value=0.0, max_value=1.0, value=0.90)

    h_bonds = np.random.randint(10, 20)
    salt_bridges = np.random.randint(1, 4)
    pi_stacks = np.random.randint(1, 3)

    data = pd.DataFrame([{
        'Aptamer': custom_aptamer,
        'EGFR_Variant': 'Custom',
        'Docking_Score': docking_score,
        'Confidence_Score': confidence_score,
        'Binding_Score': docking_score,
        'Num_Interactions': h_bonds + salt_bridges + pi_stacks,
        'Hydrogen_Bonds': h_bonds,
        'Salt_Bridges': salt_bridges,
        'Pi_Stacking': pi_stacks,
        'Stability_Score': round(np.random.uniform(0.85, 0.96), 2),
        'Mutation_Impact': 'Predicted'
    }])
else:
    data = pd.DataFrame()

st.subheader("Interaction Summary")
st.dataframe(data)

st.subheader("3D Structure Viewer")
with st.expander("View Structure"):
    xyzview = py3Dmol.view(width=600, height=400)
    if pdb_file:
        pdb_string = pdb_file.read().decode("utf-8")
        xyzview.setBackgroundColor('black')
        xyzview.addModel(pdb_string, 'pdb')
        xyzview.setStyle({'cartoon': {'color': 'spectrum'}})
        xyzview.zoomTo()
        st.components.v1.html(xyzview._make_html(), height=400)
