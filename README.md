# Cocultures2018

The full analysis of the submitted paper ***Coculturing bacteria leads to reduced phenotypic heterogeneities*** by Jasmine Heyse, Benjamin Buysschaert, Ruben Props, Peter Rubbens, Andre Skirtach, Willem Waegeman and Nico Boon.

The analysis pipeline can be found in _AnalysisCocultures.html_. Before starting the analysis please unzip Ramanfiles.zip and store them in _Data\Ramanfiles_. Download the FCM data from FlowRepository (accession ID FR-FCM-ZYWN) and store them in a folder named _FCSfiles_. Metadata should be stored in separate files in the _Data_ folder.

The final file structure should be: 

```
├── AnalysisCocultures.Rmd
├── AnalysisCocultures.html
├── AnalysisCocultures.md
├── bibliography.bib
├── CitationStyle.csl
├── /Data
	├── Metadata_FCM.csv
    ├── Metadata_Raman.csv
    ├── /FCSfiles
    └── /Ramanfiles
```