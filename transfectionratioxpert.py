import json
import os
from datetime import datetime

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st
from scipy.optimize import minimize


culture_vessel_options_lipo = {"96-well": 25, "24-well": 50, "12-well": 100, "6-well/35-mm": 250,
                               "60-mm/flask 25 cm2": 500, "10-cm/flask 75 cm2": 1500}
culture_vessel_options_jetprime = {"96-well": 5, "24-well": 50, "12-well": 75, "6-well/35-mm": 200,
                                   "60-mm/flask 25 cm2": 200, "10-cm/flask 75 cm2": 500, "15-cm/flask 175cm2": 1000}
calcul = []

# Settings for Streamlit page
st.set_page_config(
    page_title="TRXpert",
    page_icon="🧬",
    layout="wide")

# Main page
st.sidebar.title('🧬 Transfection Ratio Xpert')

df = pd.DataFrame(
    [
        {"Plasmid/Vector/RNA": "pCDNA3.1", "µg/µL": 1.5},
        {"Plasmid/Vector/RNA": "pSBBi-GP", "µg/µL": 2}
    ]
)
df = st.sidebar.data_editor(df, key="vector", num_rows="dynamic")
vector_table = df.to_dict(orient='records')

num_columns = st.sidebar.number_input('Mix', min_value=1, max_value=None, value=1, step=1, help=None, key="mix")
columns = st.columns(num_columns)

for i in range(1, num_columns + 1):
    with columns[i % num_columns - 1]:
        vector_per_mix = st.number_input('Vector per mix', min_value=1, max_value=len(df), value=1, step=1, help=None,
                                         key=f"vector_per_mix{i}")
        amount_of_dna = st.number_input('Amount of total DNA (µg)', min_value=0.0, max_value=10.0, value=1.0, step=0.1,
                                        help=None, key=f'amount_dna{i}')
        transfection_type = st.radio("Transfection type", ["Lipofectamine (2000/3000)", 'jetPRIME'], key=f"transfection_type{i}", horizontal=True)
        dna_ratio = st.number_input("DNA ratio", min_value=0.1, max_value=None, value=1.0, step=0.1, help=None,
                                    key=f"dna_ratio{i}")
        reagent_ratio = st.number_input("Reagent ratio", min_value=0.1, max_value=None, value=2.0, step=0.1, help=None,
                                        key=f"reagent_ratio{i}")
        ratio = dna_ratio / reagent_ratio
        culture_vessel_value = culture_vessel_options_lipo[
            st.selectbox("Culture Vessel", list(culture_vessel_options_lipo.keys()), key=f"culture_vessel{i}")] if transfection_type == "Lipofectamine (2000/3000)"\
            else culture_vessel_options_jetprime[st.selectbox("Culture Vessel", list(culture_vessel_options_jetprime.keys()), key=f"culture_vessel{i}")]
        number_wells = st.number_input("Number of well(s)", min_value=1, max_value=None, value=1, step=1, help=None, key=f"number_well{i}")
        for j in range(1, vector_per_mix + 1):
            vector_selected = st.selectbox(f'Vector {j}', df["Plasmid/Vector/RNA"], key=f'vector_selected{i}-{j}')
            sum_amount_dna_selected_same_i = sum(amount_dna_selected for item in calcul if item[0] == i)
            max_slider = amount_of_dna - sum_amount_dna_selected_same_i
            if max_slider > 0:
                amount_dna_selected = st.slider(f'Amount of {vector_selected}', 0.0, max_slider, step=0.1,
                                                key=f'amount_vector_selected{i}-{j}')
            else:
                st.warning(f"Amount of {vector_selected} not available. The previous plasmid(s)/vector(s)/RNA already uses all the DNA required for transfection.")
                amount_dna_selected = 0
            calcul.append([i, amount_of_dna, transfection_type, ratio, culture_vessel_value, number_wells, vector_selected, amount_dna_selected])


results_dict = {}
for i, amount_of_dna, transfection_type, ratio, culture_vessel_value, number_wells, vector_selected, amount_dna_selected in calcul:
    if vector_selected not in results_dict:
        results_dict[vector_selected] = {
            "Plasmid/Vector/RNA": vector_selected,
            f"Mix {i} (µL)": (amount_dna_selected / next(item['µg/µL'] for item in vector_table if item["Plasmid/Vector/RNA"] == vector_selected)) * number_wells
        }
    else:
        results_dict[vector_selected][f"Mix {i} (µL)"] = (amount_dna_selected / next(item['µg/µL'] for item in vector_table if item["Plasmid/Vector/RNA"] == vector_selected)) * number_wells

    if transfection_type not in results_dict:
        results_dict[transfection_type] = {
                    "Plasmid/Vector/RNA": transfection_type, f"Mix {i} (µL)": (amount_of_dna / ratio) * number_wells}
    else:
        results_dict[transfection_type][f"Mix {i} (µL)"] = (amount_of_dna / ratio) * number_wells

    if transfection_type == "Lipofectamine (2000/3000)":
        if "Culture Vessel Lipo 1" not in results_dict:
            results_dict["Culture Vessel Lipo 1"] = {
                "Plasmid/Vector/RNA": "OptiMEM w/ Plasmid/Vector/RNA", f"Mix {i} (µL)": culture_vessel_value * number_wells}
            results_dict["Culture Vessel Lipo 2"] = {
                "Plasmid/Vector/RNA": "OptiMEM w/ Lipofectamine", f"Mix {i} (µL)": culture_vessel_value * number_wells}
        else:
            results_dict["Culture Vessel Lipo 1"][f"Mix {i} (µL)"] = culture_vessel_value * number_wells
            results_dict["Culture Vessel Lipo 2"][f"Mix {i} (µL)"] = culture_vessel_value * number_wells

    elif transfection_type == "jetPRIME":
        if "Culture Vessel JP" not in results_dict:
            results_dict["Culture Vessel JP"] = {"Plasmid/Vector/RNA": "jetPRIME Buffer", f"Mix {i} (µL)": culture_vessel_value * number_wells}
        else:
            results_dict["Culture Vessel JP"][f"Mix {i} (µL)"] = culture_vessel_value * number_wells

for i, amount_of_dna, transfection_type, ratio, culture_vessel_value, number_wells, vector_selected, amount_dna_selected in calcul:
    if transfection_type == "Lipofectamine (2000/3000)":
        if "Culture Vessel Lipo 1" in results_dict:
            culture_entry = results_dict.pop("Culture Vessel Lipo 1")
            results_dict["Culture Vessel Lipo 1"] = culture_entry

        if "Culture Vessel Lipo 2" in results_dict:
            culture_lipo_entry = results_dict.pop("Culture Vessel Lipo 2")
            results_dict["Culture Vessel Lipo 2"] = culture_lipo_entry

        if transfection_type in results_dict:
            reagent_entry = results_dict.pop(transfection_type)
            results_dict[transfection_type] = reagent_entry

    else:
        if "Culture Vessel JP" in results_dict:
            culture_entry = results_dict.pop("Culture Vessel JP")
            results_dict["Culture Vessel JP"] = culture_entry

        if transfection_type in results_dict:
            reagent_entry = results_dict.pop(transfection_type)
            results_dict[transfection_type] = reagent_entry


df_results = pd.DataFrame(list(results_dict.values()))

st.dataframe(df_results, hide_index=True)