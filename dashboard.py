import streamlit as st
import py3Dmol
import pandas as pd

st.set_page_config(layout="wide")
st.title("EGFR-Aptamer Interaction Dashboard")

# 3D Viewer Section
st.subheader("1. EGFR Protein Structure")
pdb_id = "1M17"  # You can change this later
view = py3Dmol.view(query=f"pdb:{pdb_id}")
view.setStyle({'cartoon': {'color': 'spectrum'}})
view.zoomTo()
view.show()

# Aptamer Upload
st.subheader("2. Upload Aptamer Sequence (FASTA or TXT)")
aptamer_file = st.file_uploader("Upload Aptamer", type=["fasta", "txt"])

if aptamer_file:
    content = aptamer_file.read().decode("utf-8")
    
    if content.startswith(">"):
        lines = content.strip().split("\n")
        header = lines[0].replace(">", "")
        sequence = "".join(lines[1:])
    else:
        header = "Custom Aptamer"
        sequence = content.strip()

    # Validate sequence (basic ACGT/N check)
    import re
    if not re.match("^[ACGTUNacgtun]+$", sequence):
        st.warning("This doesn't look like a valid DNA/RNA sequence.")
    else:
        st.markdown(f"**Aptamer Name:** `{header}`")
        st.text_area("Aptamer Sequence", sequence, height=150)
        st.success(f"Aptamer '{header}' uploaded successfully!")

        # Store for later use
        st.session_state["aptamer_seq"] = sequence
        st.session_state["aptamer_name"] = header
# Simulate Docking & Interaction
st.subheader("3. Simulated Docking & Interaction")
if st.button("Run Docking (Simulated)"):
    st.success("Docking Complete. Showing mock interaction results.")
    df = pd.read_csv("mock_results/aptamer_interactions.csv")
    st.dataframe(df)

    best = df.sort_values("Num_Bonds", ascending=False).iloc[0]
    st.info(f"Top Aptamer: {best['Aptamer']} with {best['Num_Bonds']} bonds.")
