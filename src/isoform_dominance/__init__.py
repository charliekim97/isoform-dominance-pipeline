"""Is this functional isoform-class comparison measurable by short reads? Group-level
identifiability, quantification and testing from bulk RNA-seq

Modules
-------
annotate        gene symbol -> proposed isoform groups (Ensembl)
identifiability are two isoform groups distinguishable by short reads?
extract         Salmon quant.sf -> per-donor isoform-group TPM
stats           paired Wilcoxon + figure
contamination   marker-based contamination control
"""
__version__ = "2.9.0.dev0"
