import streamlit as st
import pandas as pd
import py3Dmol
import matplotlib.pyplot as plt
import base64
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from io import StringIO

# Sample training data (mocked for ML model)
train_data = pd.DataFrame({
    'Hydrogen_Bonds': [10, 15, 12, 14, 13],
    'Salt_Bridges': [2, 3, 2, 2, 3],
    'Pi_Stacking': [1, 2, 2, 2, 1],
    'Binding_Score': [-180, -215.6, -190, -205.2, -198.5]
})
X_train = train_data[['Hydrogen_Bonds', 'Salt_Bridges', 'Pi_Stacking']]
y_train = train_data['Binding_Score']
model = RandomForestRegressor()
model.fit(X_train, y_train)

st.set_page_config(layout="wide")
st.title("RNA Aptamer Docking & Analysis Resource - R.A.D.A.R.")

# Upload section
st.sidebar.header("Upload Custom EGFR Variant")
pdb_file = st.sidebar.file_uploader("Upload PDB File", type=['pdb'])
custom_aptamer = st.sidebar.text_input("Enter Custom Aptamer Label", "APT-CUSTOM")

# Ensure session state holds data even before any upload
if 'uploaded_data' not in st.session_state:
    st.session_state.uploaded_data = pd.DataFrame(columns=[
        'Aptamer', 'EGFR_Variant', 'Docking_Score', 'Confidence_Score',
        'Binding_Score', 'Num_Interactions', 'Hydrogen_Bonds',
        'Salt_Bridges', 'Pi_Stacking', 'Stability_Score', 'Mutation_Impact'
    ])

data = st.session_state.uploaded_data.copy()

# Sidebar filters
aptamer_list = data['Aptamer'].unique().tolist()
selected_aptamers = st.sidebar.multiselect("Select Aptamer(s)", aptamer_list, default=aptamer_list)
variants = st.sidebar.multiselect("Select EGFR Variant(s)", ['Wild-Type', 'L858R', 'T790M', 'Custom'], default=['Wild-Type', 'L858R', 'T790M'])
score_threshold = st.sidebar.slider("Minimum Binding Score (stronger = more negative)", -250, -150, -220)

if pdb_file:
    # Simulate interaction analysis
    h_bonds = np.random.randint(10, 20)
    salt_bridges = np.random.randint(1, 4)
    pi_stacks = np.random.randint(1, 3)

    docking_score = st.number_input("Enter docking score (e.g., from HDOCK)", value=-200.0)
    confidence_score = st.number_input("Enter confidence score (0–1)", min_value=0.0, max_value=1.0, value=0.90)

    new_row = pd.DataFrame([{
        'Aptamer': custom_aptamer,
        'EGFR_Variant': 'Custom',
        'Docking_Score': docking_score,
        'Confidence_Score': confidence_score,
        'Binding_Score': docking_score,  # Or separate if needed
        'Num_Interactions': h_bonds + salt_bridges + pi_stacks,
        'Hydrogen_Bonds': h_bonds,
        'Salt_Bridges': salt_bridges,
        'Pi_Stacking': pi_stacks,
        'Stability_Score': round(np.random.uniform(0.85, 0.96), 2),
        'Mutation_Impact': 'Predicted'
    }])

    # Save or append
    if 'uploaded_data' not in st.session_state:
        st.session_state.uploaded_data = pd.DataFrame(columns=[
        'Aptamer', 'EGFR_Variant', 'Docking_Score', 'Confidence_Score',
        'Binding_Score', 'Num_Interactions', 'Hydrogen_Bonds',
        'Salt_Bridges', 'Pi_Stacking', 'Stability_Score', 'Mutation_Impact'
    ])

data = st.session_state.uploaded_data.copy()

# Filtered data
filtered_data = data[
    (data['Aptamer'].isin(selected_aptamers)) &
    (data['EGFR_Variant'].isin(variants)) &
    (data['Binding_Score'] <= score_threshold)
]

# Display table
st.subheader("Interaction Summary")
st.write(data.columns)  # 🔍 DEBUG: Show available columns
st.dataframe(data[['Aptamer', 'EGFR_Variant', 'Docking_Score', 'Confidence_Score', 'Binding_Score', 'Stability_Score']])

# 3D Structure Viewer
st.subheader("3D Structure Viewer")
with st.expander("View Structure"):
  xyzview = py3Dmol.view(width=400, height=400)
xyzview.setBackgroundColor('white')
xyzview.addModel(pdb_string, 'pdb')


# Heatmap
st.subheader("Binding Score Heatmap")
heatmap_data = filtered_data.pivot(index='EGFR_Variant', columns='Aptamer', values='Docking_Score')
if not heatmap_data.empty:
    fig, ax = plt.subplots()
    cax = ax.matshow(heatmap_data, cmap='coolwarm')
    plt.xticks(range(len(heatmap_data.columns)), heatmap_data.columns)
    plt.yticks(range(len(heatmap_data.index)), heatmap_data.index)
    plt.colorbar(cax)
    st.pyplot(fig)
else:
    st.warning("No data to display for the current filters.")

# Base style for all residues
xyzview.setStyle({'cartoon': {'color': 'spectrum'}})

xyzview.setStyle({'cartoon': {'color': 'white'}})  # base for everything
xyzview.addStyle({'chain': 'A'}, {'cartoon': {'color': 'cyan'}})
xyzview.addStyle({'chain': 'B'}, {'cartoon': {'color': 'orange'}})
xyzview.addStyle({'resn': 'ASN'}, {'cartoon': {'color': 'blue'}})
xyzview.addStyle({'resn': 'ALA'}, {'cartoon': {'color': 'red'}})
xyzview.addStyle({'resn': 'LIG'}, {'stick': {'color': 'magenta'}})
xyzview.zoomTo()
st.components.v1.html(xyzview._make_html(), height=400)


st.write("Boxplot columns:", filtered_data.columns.tolist())
# Boxplot of Docking Scores by EGFR Variant
st.subheader("Docking Score Distribution by EGFR Variant")
if not filtered_data.empty:
    fig_box, ax_box = plt.subplots()
    filtered_data.boxplot(column='Docking_Score', by='EGFR_Variant', ax=ax_box)
    plt.suptitle("")  # Clean up automatic title
    ax_box.set_title("Boxplot of Docking Scores")
    ax_box.set_ylabel("Docking Score")
    st.pyplot(fig_box)
else:
    st.info("No data available to show boxplot.")
# Bar Chart of Interactions
st.subheader("Interaction Breakdown")
if not filtered_data.empty:
    fig2, ax2 = plt.subplots()
    bond_types = ['Hydrogen_Bonds', 'Salt_Bridges', 'Pi_Stacking']
    interaction_counts = filtered_data[bond_types].sum()
    interaction_counts.plot(kind='bar', ax=ax2)
    plt.ylabel("Total Count")
    st.pyplot(fig2)

# Mutation Sensitivity
st.subheader("Mutation Sensitivity Report")
if not filtered_data.empty:
    fig3, ax3 = plt.subplots()
    filtered_data.groupby('EGFR_Variant')['Stability_Score'].mean().plot(kind='bar', ax=ax3)
    plt.ylabel("Avg. Stability Score")
    plt.title("Stability Across Variants")
    st.pyplot(fig3)

    st.markdown("**Mutation Impact Summary**")
    impact_table = filtered_data[['EGFR_Variant', 'Mutation_Impact', 'Stability_Score']]
    st.dataframe(impact_table)
#scatter plot
st.subheader("Docking Score vs Stability")
if not filtered_data.empty:
    fig_scatter, ax_scatter = plt.subplots()
    ax_scatter.scatter(filtered_data['Docking_Score'], filtered_data['Stability_Score'], color='blue')
    ax_scatter.set_xlabel("Docking Score")
    ax_scatter.set_ylabel("Stability Score")
    ax_scatter.set_title("Docking Score vs Stability")
    st.pyplot(fig_scatter)
else:
    st.info("No data to show scatter plot.")

# Data Export
st.subheader("Export Data")
@st.cache_data
def convert_df(df):
    return df.to_csv(index=False).encode('utf-8')

csv = convert_df(filtered_data)
st.download_button(
    label="Download Filtered Data as CSV",
    data=csv,
    file_name='aptamer_egfr_data.csv',
    mime='text/csv'
)

st.markdown("---")
st.markdown("Developed as part of an advanced computational biology project on aptamer-EGFR targeting in glioblastoma.")
