# chemical_conversion
## Small toolbox to retrieve InChIKeys from various chemical representations 
Collection of Python scripts using the ChEMBL structure curation pipeline and RDKit. 
For further details check the GitHub repo of the ChEMBL structure curation pipeline (https://github.com/chembl/ChEMBL_Structure_Pipeline/wiki) or the original publication [[1]](#1). 

## Motivation
Chemicals, in particular small molecule drugs used in academic research as provided by various vendors, are often named ambiguously. Furthermore, different chemicals with the same parental structure (e.g. different salts of the same molecule) are annotated as different chemicals in commercial libraries - often with varying meta annotations (e.g. drug targets). The described issues are particularly severe if chemicals or annotations are retrieved from different sources. Thus, an ambiguous identification of chemicals is essential to harmonize libraries and annotations. We here use a published pipeline to get the standardized parental structure of a chemical and then retrieve the InChIKey representation. InChIKey can then be used as robust identifiers for chemicals.

## References
<a id="1">[1]</a> 
Bento, A.P., Hersey, A., Félix, E. et al. An open source chemical structure curation pipeline using RDKit. J Cheminform 12, 51 (2020). https://doi.org/10.1186/s13321-020-00456-1
