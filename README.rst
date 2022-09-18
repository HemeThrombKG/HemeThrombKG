HemeThrombKG
============

Modeling the link between labile heme and thrombosis.
*******************************************************

This repository contains data, notebooks and scripts to perform the analyses and generate figures as presented in the
paper, *Exploring the complex network of heme-triggered effects on the blood coagulation system.*

`HemeThrombKG <https://github.com/sarahbeenie/hemeThrombKG/blob/main/bel_files/hemeThrombKG.bel>`_ is a knowledge graph
which models the interactions between heme and its interference in coagulation processes, as summarized in a review by
Hopp *et al.* (2021) [1]_.


Folder Structure
----------------

The folder structure of this repository is organized as follows:

::

    HemeThrombKG
    ├── LICENSE
    ├── README.rst
    ├── analyses  # Jupyter notebooks and script for analyses and visualizations
    │   ├── README.md
    │   ├── cellular_analysis.ipynb
    │   ├── extracellular_protein_overlay.ipynb
    │   ├── intracellular_protein_overlay.ipynb
    │   ├── network_utils.py
    │   ├── node_frequency_analysis.ipynb
    │   ├── node_frequency_overlay_extracellular.ipynb
    │   ├── node_frequency_overlay_intracellular.ipynb
    │   ├── protein_analysis_extracellular.ipynb
    │   └── protein_analysis_intracellular.ipynb
    ├── bel_files  # BEL files to generate HemeThrombKG and pathway networks
    │   ├── coagulation_pathways.bel
    │   ├── common_pathway_reactome.bel
    │   ├── extrinsic_pathway_reactome.bel
    │   ├── hemeThrombKG.bel
    │   ├── intrinsic_pathway_reactome.bel
    │   ├── platelet_activation_kegg.bel
    │   └── plug_formation_reactome.bel
    ├── figures  # Paper figures
    │   ├── extracellular proteins
    │   ├── intracellular proteins
    │   └── node frequency
    └── node_type_files  # Node label files
        ├── cell_terms.tsv
        ├── extracellular_proteins.tsv
        └── intracellular_proteins.tsv



References
----------

.. [1] Mubeen S., Domingo-Fernandez D., Díaz del Ser S., Solanki D., Kodamullil A.T., Hofmann-Apitius M.,
    Hopp M.T., and Imhof D. [*]_ (2022). Exploring the complex
    network of heme-triggered effects on the blood coagulation system. *Preprint.*

.. [2] Hopp M.T., and Imhof D. (2021). Disease networks. `Linking Labile Heme with Thrombosis <https://doi.org/10.3390/jcm10030427>`_. Journal of Clinical Medicine. 10(3), 427.


.. [*] Hopp M.T., and Imhof D. have contributed equally to this work.

