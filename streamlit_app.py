import streamlit as st

# from streamlit import streamlit_js_eval
# import streamlit-js-eval
import pandas as pd
import altair as alt
import re
import numpy as np
from sklearn.linear_model import LinearRegression


st.set_page_config(page_title="CRISPRi_vis", page_icon=":test_tube:", layout="wide")


@st.cache_data
def load_data(data_url):
    df = pd.read_csv(data_url)
    df.index = df["ID"]
    return df

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

data_url = "https://github.com/aphren/Hren-2024-scripts/raw/main/library_data_with_loci.csv"
data = load_data(data_url)

## extract a dictionary to retrieve gene descriptions
gene_description_unique = data.drop_duplicates(subset="gene", keep="last")
gene_description_dict = pd.Series(
    gene_description_unique["product"].values, index=gene_description_unique["gene"]
).to_dict()

# growth_condition = '22R'

st.sidebar.markdown(
    """<center>
    <h2>Welcome to our CRISPRi data visualization tool! <br> <a href="https://www.biorxiv.org/content/10.1101/2024.05.20.595006v1">Hren et al. 2024</a> <h2>
    Enter individual genes on the main page, or observe the behavior of select protein groups
    <br>
""",
    unsafe_allow_html=True,
)

protein_complexes = [   
    "NDH-1M (core subunits)",
    "Photosystem II",
    "Photosystem I",
    "Phycobilisome",
    "tRNA synthetase",
    "Purine metabolism",
    "Circadian rhythm",
    "Glycolysis",
    "Ribosome",
    "Transcription",
]

all_genes = data["gene"].unique().tolist()
## Select the protein complex or functional group of interest
complex_selection = st.sidebar.selectbox(
    "Select protein complex or gene group:",
    protein_complexes,
    placeholder="Choose from dropdown",
    index=None,
)

### Create a list so people can look for gene names
gene_dropdown = st.sidebar.selectbox(
    "Gene descriptions:", all_genes, placeholder="Type or select", index=None
)
if gene_dropdown:
    gene_description = gene_description_dict.get(gene_dropdown, [])
    st.sidebar.write(gene_description)


color_scheme_box = st.sidebar.checkbox("Change color scheme")
if color_scheme_box:
    color_scheme = "category10"
else:
    color_scheme = "tableau20"


ndh1m = [
    "ndhA",
    "ndhB",
    "ndhC",
    "ndhE",
    "ndhG",
    "ndhH",
    "ndhI",
    "ndhJ",
    "ndhK",
    "ndhL",
    "ndhM",
    "ndhN",
    "ndhO",
    "ndhP",
]
PS2 = [
    "psbA",
    "psbA-II",
    "psbB",
    "psbC",
    "psbD",
    "psbE",
    "psbF",
    "psbH",
    "psbJ",
    "psbK",
    "psbL",
    "psbM",
    "psbN",
    "psbO",
    "psbT",
    "psbU",
    "psbV",
    "psbW",
    "psbW2",
    "psbY",
    "psbZ",
]
PS1 = [
    "psaA",
    "psaB",
    "psaC",
    "psaD",
    "psaE",
    "psaF",
    "psaI",
    "psaJ",
    "psaK",
    "psaL",
    "psaM",
]
PBS = [
    "cpcA",
    "cpcB",
    "cpcC",
    "cpcD",
    "cpcE",
    "cpcF",
    "cpcG",
    "apcA",
    "apcB",
    "apcC",
    "apcD",
    "apcE",
    "apcF",
]
synthetases = [
    "alaS",
    "asnS",
    "aspS",
    "gltX",
    "glyQ",
    "hisS",
    "ileS",
    "lysS",
    "gatB",
    "pheS",
    "proS",
    "trpS",
    "serS",
    "thrS",
    "tyrS",
    "valS",
    "argS",
    "cysS",
    "leuS",
    "metG",
]
purines = [
    "purA",
    "purB",
    "purC",
    "purD",
    "purE",
    "purF",
    "purH",
    "purK",
    "purL",
    "purM",
    "purN",
    "purQ",
    "purS",
    "purT",
    "purU",
]
clock = ["kaiA", "kaiB", "kaiC", "rpaA", "sasA"]
glycolysis = ["gap", "pgk", "tpiA", "fba", "pfkA", "pgi", "gpm", "eno", "glk", "pyk"]
ribosomes = ["rplA", "rplB", "rplC", "rplL", "rpmA", "rpmB", "rpmC", "rpmG"]
transcription = [
    "rpoA",
    "rpoB",
    "rpoC1",
    "rpoC2",
    "rpoD",
    "rpoZ",
    "gyrA",
    "gyrB",
    "pcrA",
]

complex_dict = {
    "NDH-1M (core subunits)": ndh1m,
    "Photosystem II": PS2,
    "Photosystem I": PS1,
    "Phycobilisome": PBS,
    "tRNA synthetase": synthetases,
    "Purine metabolism": purines,
    "Circadian rhythm": clock,
    "Glycolysis": glycolysis,
    "Ribosome": ribosomes,
    "Transcription": transcription,
}


if complex_selection:
    gene_list = complex_dict.get(complex_selection, [])
else:
    gene_input = st.text_input(
        "Enter gene names: separate multiple genes with a space", key="gene_user_input"
    )
    gene_list = gene_input.split()

condition_list = [
    "22°C, Red",
    "22°C, Blue",
    "22°C, White",
    "37°C, Red",
    "37°C, Blue",
    "37°C, Diurnal",
]

condition_dict = {
    "22°C, Red": "22R",
    "22°C, Blue": "22B",
    "22°C, White": "22W",
    "37°C, Red": "37R",
    "37°C, Blue": "37B",
    "37°C, Diurnal": "37di",
}

growth_condition = st.selectbox("**Select growth condition:**", condition_list)


if gene_list:
    chart_data = pd.DataFrame()
    ## Define relevant column
    condition_column = condition_dict[growth_condition] + "_LFC"

    ## Extract condition info for axis titles
    regex_pattern = r"(\d{2})([A-Za-z]{1,2})"
    match = re.match(regex_pattern, condition_dict[growth_condition])
    colors = {"B": "Blue", "R": "Red", "W": "White", "di": "Diurnal"}
    if match:
        temperature = match.group(1)
        light = colors[match.group(2)]

    condition_title = f"Log₂ Fold Change ({temperature}°C, {light})"

    for gene in gene_list:
        temp_chart_data = data[data["gene"] == gene][
            ["Spacer", "37W_LFC", condition_column]
        ]
        temp_chart_data["gene"] = gene

        # fit linear regression to each displayed gene
        try:
            #model = LinearRegression()
            X = temp_chart_data[["37W_LFC"]].values.reshape(-1, 1)
            y = temp_chart_data[condition_column].values
            temp_chart_data["regression_y"] = moving_average(y,2)
            
            #model.fit(X, y)
            #temp_chart_data["regression_y"] = model.predict(X)
        except:
            pass
        chart_data = pd.concat([chart_data, temp_chart_data])

    legend_selection = alt.selection_point(fields=["gene"], bind="legend")

    base_chart = (
        alt.Chart(chart_data)
        .mark_point()
        .encode(
            x=alt.X(
                "37W_LFC:Q", title="Log₂ Fold Change (37°C, White)"
            ),  # Replace with your actual column name for x-axis
            y=alt.Y(
                condition_column, type="quantitative", title=condition_title
            ),  # Replace with your actual column name for y-axis
            color=alt.Color(
                "gene:N", scale=alt.Scale(scheme=color_scheme)
            ),  # Darker color scheme, https://vega.github.io/vega/docs/schemes/
            opacity=alt.condition(legend_selection, alt.value(1), alt.value(0.2)),
        )
        .properties(
            height=550,
            title={
                "text": ["Gene Fold Change Comparison"],
                "subtitle": ["Click gene name in legend to highlight"],
            },
        )
        .add_params(legend_selection)
    )
    # Add regression lines
    regression_lines = (
        alt.Chart(chart_data)
        .mark_line()
        .encode(
            x=alt.X("37W_LFC:Q"),
            y=alt.Y("regression_y:Q"),  # Adjust as needed
            color=alt.Color("gene:N", scale=alt.Scale(scheme=color_scheme)),
            tooltip=["gene", "regression_y:Q"],
        )
    )

    zero_line_y = (
        alt.Chart(pd.DataFrame({"y": [0]}))
        .mark_rule(color="red", strokeDash=[3, 3])
        .encode(y="y:Q")
    )
    zero_line_x = (
        alt.Chart(pd.DataFrame({"x": [0]}))
        .mark_rule(color="red", strokeDash=[3, 3])
        .encode(x="x:Q")
    )

    linear_regressions = st.sidebar.checkbox("Add linear regressions")

    if linear_regressions:
        chart = regression_lines + base_chart + zero_line_y + zero_line_x
    else:
        chart = base_chart + zero_line_y + zero_line_x

    col1, col2 = st.columns([3, 2], gap="large")
    with col1:
        st.altair_chart(chart, use_container_width=True)

    with col2:
        st.write("**Guide data:**")
        st.write(chart_data)
else:
    st.write("Please enter gene names or select a protein complex from the sidebar.")
