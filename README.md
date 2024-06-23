# ModuleInference
Python code to reproduce the correlation-based inference of regulatory signaling and transcriptional modules from time-series molecular data, as published in:

Flower CT*, Liu C*, Chuang H-Y, Ye X, Cheng H, Heath JR**, Wei W**, White FM**. Signaling and transcriptional dynamics underlying early adaptation to oncogenic BRAF inhibition. _bioRxiv_. 2024;581004.

Run module_inference.py in the same directory as the published phosphoproteomic and transcriptomic data (**Supplementary_Table_S1.xlsx** & **Supplementary_Table_S2.xlsx**) and the **gene_set_libraries** directory. This will generate a number of intermediate csv-formatted files, which are used by module_integration.py to graph the dynamics of each regulatory module and to calculate the pairwise meta-correlation between modules.

If you make use of this code or framework for analysis of time-series molecular data, please cite our corresponding manuscript above.
