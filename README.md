# DIAMOnD-Based Disease Module Identification for Alzheimer's Disease

## Project overview

This project applies DIAMOnD (Disease Module Detection) to identify proteins that may be related to Alzheimer's disease through their positions in a human protein-protein interaction (PPI) network.

The analysis starts with known Alzheimer's disease-associated seed genes and ranks additional proteins using the DIAMOnD algorithm. It then explores the effect of the DIAMOnD alpha parameter and adds complementary analyses:
- Louvain community detection on the PPI network
- Comparison of predictions generated with alpha = 1, 2 and 10
- Comparison with Alzheimer's disease target-association evidence from Open Targets
- Pathway enrichment of the alpha = 10 candidate genes using Enrichr

## Research question

Can DIAMOnD identify additional candidate proteins associated with Alzheimer's disease from a PPI network, and how does changing the alpha parameter affect the predicted candidates and their supporting biological evidence?

## Data sources

The project uses the supplementary datasets associated with the DIAMOnD research paper:

- Human interactome: https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1004120.s003&type=supplementary
- Disease-gene associations: https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1004120.s004&type=supplementary
- Original DIAMOnD implementation: https://github.com/barabasilab/DIAMOnD

The notebook extracts the Alzheimer's disease seed genes from the disease-gene supplementary file and retains seed genes present in the interactome. The resulting prepared inputs are `PPI.txt` and `seed_genes.txt`.

External evidence sources used in the notebook:
- MyGene.info for mapping gene symbols to Ensembl identifiers: https://mygene.info/
- Open Targets Platform for Alzheimer's disease target associations: https://platform.opentargets.org/
- Enrichr for pathway enrichment: https://maayanlab.cloud/Enrichr/

## Repository structure

```text
DIAMOnD-Alzheimers-Project/
├── README.md
├── requirements.txt
├── download_data.py
├── code/
│   ├── DIAMOnD.py
│   └── DIAMOnD_py3.py
├── notebooks/
│   └── SE23UCAB010_DIAMOnD_Project.ipynb
├── results/
│   ├── alzheimer_diamond_alpha1.txt
│   ├── alzheimer_diamond_alpha2.txt
│   ├── alzheimer_diamond_alpha10.txt
│   ├── DIAMOnD_Alpha10_All_Pathways.csv
│   └── DIAMOnD_Alpha10_Significant_Pathways.csv
└── figures/
    ├── Alzheimer_DIAMOnD_OpenTargets_Top20.png
    └── Alzheimer_Alpha_Validation_Comparison.png
```

The `download_data.py` script downloads the original supplementary datasets and creates the prepared `PPI.txt` and `seed_genes.txt` input files. Raw data files do not need to be committed to GitHub.

## How to run

### 1. Install dependencies

Use Python 3. Run:

```bash
pip install -r requirements.txt
```

### 2. Download and prepare the input data

From the repository root, run:

```bash
python download_data.py
```

This creates:
- `data/raw/original_interactome.tsv`
- `data/raw/original_disease_genes.tsv`
- `data/processed/PPI.txt`
- `data/processed/seed_genes.txt`

### 3. Run DIAMOnD

The notebook commands used during the analysis are equivalent to:

```bash
python code/DIAMOnD_py3.py data/processed/PPI.txt data/processed/seed_genes.txt 100 1 results/alzheimer_diamond_alpha1.txt
python code/DIAMOnD_py3.py data/processed/PPI.txt data/processed/seed_genes.txt 100 2 results/alzheimer_diamond_alpha2.txt
python code/DIAMOnD_py3.py data/processed/PPI.txt data/processed/seed_genes.txt 100 10 results/alzheimer_diamond_alpha10.txt
```

The exact input/output paths may need to be adjusted if you run the notebook in a different directory. Open the notebook in Jupyter or Google Colab and execute its cells in order for the full downstream analysis, including external API queries.

## Analysis and results

The notebook:
1. Loads the human PPI network and Alzheimer's disease seed genes.
2. Runs DIAMOnD to obtain 100 candidate proteins for alpha = 1, 2 and 10.
3. Compares candidate overlap and Jaccard similarity across alpha settings.
4. Uses Louvain community detection to examine the distribution of seed genes and candidate proteins in network communities.
5. Maps candidate gene symbols to Ensembl identifiers with MyGene.info.
6. Compares candidates with Open Targets Alzheimer's disease associations.
7. Submits the alpha = 10 candidate list to Enrichr and retrieves KEGG 2021 Human and Reactome 2022 enrichment results.
8. Produces figures summarizing Open Targets scores and average validation scores across alpha settings.

### Reported findings from the completed analysis

- The PPI network contained 141,296 interactions and approximately 13,460 proteins.
- The analysis used 29 Alzheimer's disease seed genes.
- DIAMOnD generated 100 predictions for each of alpha = 1, 2 and 10.
- Candidate overlap changed with alpha: the reported Jaccard similarities were 48.15% for alpha 1 vs alpha 2, 19.05% for alpha 1 vs alpha 10, and 33.33% for alpha 2 vs alpha 10.
- Louvain analysis identified 118 communities. The largest community by seed overlap contained 14 of the 29 seed genes.
- Open Targets matching after gene-identifier handling found 93/100 alpha 1 candidates, 95/100 alpha 2 candidates and 98/100 alpha 10 candidates. The corresponding mean Open Targets scores reported in the analysis were 0.0735, 0.1276 and 0.1498.
- Enrichr analysis of the alpha = 10 candidate set returned pathway results from KEGG 2021 Human and Reactome 2022. The notebook reported 616 pathways with adjusted p-value below 0.05.

These are results from this particular dataset, parameter choices and analysis workflow. Open Targets association scores are evidence-aggregation scores, not probabilities that a gene causes Alzheimer's disease. Pathway enrichment indicates over-representation in the selected gene-set databases; it does not prove that every named pathway is causally involved in Alzheimer's disease.

## Limitations

- DIAMOnD prioritizes proteins using network topology and the supplied seed set; predictions are hypotheses, not proof of disease causality.
- Results depend on the PPI network, the selected seed genes and the alpha parameter.
- The three alpha settings tested are a limited parameter comparison, not an exhaustive search for an optimal value.
- External databases may be incomplete, updated over time or use different gene identifiers and aliases.
- Enrichment results depend on the selected candidate list and pathway-library versions.
- API-based validation and enrichment require an internet connection and may change as external databases are updated.

## Reproducibility notes

- The analysis notebook records the workflow and API queries.
- `download_data.py` retrieves the source supplementary datasets and prepares the DIAMOnD input files.
- Do not commit large raw datasets or API credentials. The raw data can be regenerated using the download script.
- For a fully reproducible run, use the same data versions and record the date of external API queries.

## Acknowledgements

This project uses the DIAMOnD method and supplementary data associated with the published DIAMOnD research work. Please cite the original DIAMOnD paper and method authors when reusing the algorithm or data.
