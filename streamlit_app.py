import streamlit as st
#from streamlit import streamlit_js_eval
#import streamlit-js-eval
import pandas as pd
import altair as alt
import re


st.set_page_config(
    page_title="CRISPRi_vis",
    page_icon=":test_tube:",
    layout="wide"
)


@st.cache_data
def load_data(data_url):
    df = pd.read_csv(data_url)
    df.index = df['ID']
    return df



data_url = "https://github.com/aphren/Hren-2024-scripts/raw/main/library_data.csv"
data = load_data(data_url)

## extract a dictionary to retrieve gene descriptions
gene_description_unique = data.drop_duplicates(subset='gene', keep='last')
gene_description_dict = (pd.Series(
                            gene_description_unique['product'].values, 
                            index=gene_description_unique['gene']
                        )   
                        .to_dict()
)

growth_condition = '22R'

st.sidebar.markdown("""<center>
    <h2>Welcome to our CRISPRi data visualization tool.<h2>
    Enter individual genes on the main page, or observe the behavior of select protein complexes
    <br>
""", unsafe_allow_html=True)

protein_complexes = ['NDH-1M (core subunits)', 'Photosystem II', 'Photosystem I']
all_genes = data['gene'].unique().tolist()
## Select the protein complex of interest
complex_selection = st.sidebar.selectbox(
    'Select protein complex:',
    protein_complexes,
    placeholder='Choose a complex',
    index=None
)

### Create a list so people can look for gene names
gene_dropdown = st.sidebar.selectbox(
    'Gene list:',
    all_genes,
    placeholder='Type or select',
    index=None
)

if gene_dropdown:
    gene_description = gene_description_dict.get(gene_dropdown, [])
    st.sidebar.write(gene_description)

ndh1m = ['ndhA', 'ndhB', 'ndhC', 'ndhE', 'ndhG', 'ndhH', 'ndhI', 'ndhJ', 'ndhK', 'ndhL', 'ndhM', 'ndhN', 'ndhO', 'ndhP']
PS2 = ['psbA','psbA-II','psbB','psbC','psbD','psbE','psbF','psbH','psbJ','psbK','psbL','psbM','psbN','psbO','psbT','psbU','psbV','psbW','psbW2','psbY','psbZ']
PS1 = ['psaA','psaB','psaC','psaD','psaE','psaF','psaI','psaJ','psaK','psaL','psaM']

complex_dict = {'NDH-1M (core subunits)': ndh1m,'Photosystem II': PS2,'Photosystem I': PS1}

if complex_selection:
    gene_list = complex_dict.get(complex_selection, [])
else:
    gene_input = st.text_input("Enter gene names: separate multiple genes with a space", key="gene_user_input")
    gene_list = gene_input.split()

condition_list = ['22R', '22B', '22W','37R', '37B', '37di']
growth_condition = st.selectbox(
    'Select growth condition:',
    condition_list
)

if gene_list:
    chart_data = pd.DataFrame()
    ## Define relevant column
    condition_column = growth_condition + '_LFC'

    ## Extract condition info for axis titles
    regex_pattern = r"(\d{2})([A-Za-z])"
    match = re.match(regex_pattern, growth_condition)
    colors = {'B':'Blue','R':'Red','W':'White','di':'Diurnal'}
    if match:
        temperature = match.group(1)
        light = colors[match.group(2)]

    condition_title = f'Log₂ Fold Change ({temperature}°C, {light})'


    for gene in gene_list:
        temp_chart_data = data[data['gene'] == gene][['Spacer', '37W_LFC', condition_column]]
        temp_chart_data['gene'] = gene
        chart_data = pd.concat([chart_data, temp_chart_data])


    legend_selection = alt.selection_point(fields=['gene'], bind='legend')

    base_chart = (alt.Chart(chart_data)
        .mark_point()
        .encode(
            x=alt.X('37W_LFC:Q',title='Log₂ Fold Change (37°C, White)'),  # Replace with your actual column name for x-axis
            y=alt.Y(condition_column, type='quantitative', title=condition_title),  # Replace with your actual column name for y-axis
            color='gene:N',
            opacity=alt.condition(legend_selection, alt.value(1), alt.value(0.2))
        )
        .properties(
            height=550,
            title={"text": ['Gene Fold Change Comparison'],
                   "subtitle": ['Click gene name in legend to highlight']
            }
         ).add_params(
            legend_selection
         )
    )
    zero_line = alt.Chart(pd.DataFrame({'y': [0]})).mark_rule(color='red', strokeDash=[3, 3]).encode(
        y='y:Q'
    )

    chart = base_chart + zero_line

    col1, col2 = st.columns([3, 2], gap="large")
    with col1:
        st.altair_chart(chart,use_container_width=True)
    with col2:
        st.write("Guide data:")
        st.write(chart_data)
else:
    st.write("Please enter gene names or select a protein complex from the sidebar.")
