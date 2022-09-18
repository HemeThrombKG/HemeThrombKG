# Analyses 

This folder contains analyses on the HemeThrombKG which focus on the following:

- Intracellular proteins
- Extracellular proteins
- Cell and cell-related terms
- Frequency of nodes in edges of KG

Analyses are conducted on the 
[HemeThrombKG](https://github.com/sarahbeenie/hemeThrombKG/blob/main/bel_files/hemeThrombKG.bel), HemeKG 2.0 (i.e., 
[HemeKG](https://github.com/hemekg/hemekg) and HemeThrombKG) and pathways from the KEGG and 
Reactome databases. 

## Notebooks

- **cellular_analysis:** Generate a subgraph of cell and cell-related terms within the HemeThrombKG.
- **node_frequency_analysis** Visualize the frequency of occurrences of nodes in the hemeThrombKG/HemeKG 2.0. The most
commonly occurring nodes are labelled.
- **protein_analysis:** Generate subgraphs of intracellular/extracellular nodes in the HemeThrombKG and/or HemeKG 2.0, 
pathways and their 
merged representations.
- **protein_overlay:** Overlay intracellular/extracellular nodes from the hemeThrombKG and/or HemeKG 2.0 and coagulation
pathways.
- **node_frequency_protein_overlay:** Overlay intracellular/extracellular nodes from the hemeThrombKG and/or HemeKG 2.0 
and coagulation pathways, where node sizes reflect their frequency in the KGs.

## Scripts 

- **network_utils.py**: Utils and visualization functions for HemeThrombKG and HemeKG 2.0.
