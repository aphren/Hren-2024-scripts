import streamlit as st
import pandas as pd
import altair as alt

@st.cache_data
def load_data(data_url):
    df = pd.read_csv(data_url)
    return df


data_url = "https://github.com/aphren/Hren-2024-scripts/raw/main/library_data.csv"


data = load_data(data_url)
#data_load_state.text("Done! (using st.cache_data)")

st.title("High-density CRISPRi screens reveal adaptive transcriptional gradients in cyanobacteria")




st.text_input("Enter gene names: separate multiple genes with a space (e.g. ndhG ndhE ndhC)", key="gene")

gene_list = st.session_state.gene.split()

chart_data = pd.DataFrame()
for gene in gene_list:

    temp_chart_data = (data[data['gene'] == gene]
                    [['37W_LFC','22B_LFC']]
    )
    temp_chart_data['gene'] = gene
    chart_data = pd.concat([chart_data,temp_chart_data])

    
#st.write(chart_data)
# Create Altair chart
base_chart = alt.Chart(chart_data).mark_line().encode(
    x='37W_LFC:Q',  # Replace with your actual column name for x-axis
    y='22B_LFC:Q',   # Replace with your actual column name for y-axis
    color='gene:N'
).properties(
    width=600,
    height=400,
    title='Gene Fold Change Comparison',
)


# Add a horizontal line at y = 0
zero_line = alt.Chart(pd.DataFrame({'y': [0]})).mark_rule(color='red', strokeDash=[3,3]).encode(
    y='y:Q'
)

# Layer the base chart with the zero line
chart = base_chart + zero_line

# Display the chart using Streamlit
st.altair_chart(chart)
