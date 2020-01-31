# GENE ANALYSIS 
[MEBIOINF 2019/2020]

## About the project
This project was developed in order to analyze certain genes using several Biopython modules. The main objective is to automatize the whole process of gathering information about genes and proteins and analyze it. 

## Prerequisites

- Biopython

``` pip install biopython ```

## Usage
The classes developed in ``` OV90.py ``` where created in order to retrieve gene information from online databases such as NCBI and analysing it using several tools available in Biopython. 

``` OV90.py ``` contains the following classes:

- ``` Search ``` : searches for useful literature and articles that may help understand the information retrieved from a certain gene.
- ``` RetrieveInfo ```: retrieves the gene file from NCBI and filters the information contained within it.
- ``` Blast ```: performs a _blastn_ program in order to find homologous sequences to the one being analized.
- ``` Align ```: performs alignments of either one or two+ sequences.
- ``` Clustal ```: performs multiple sequence alignment using clustalw.
- ``` Phylo ```: displays a phylogenetic tree.
- ``` SearchP ```: searches for protein information in protein related databases.
- ``` ProteinBlast ```: performs a _blastp_ program given a protein query.

To test them, a test script ``` Test.py ``` was created. So far, not all classes are implemented within ```Test.py```. As an example, the gene __STX4__ was used to test out the program. 

## Contact

- Sofia de Beir - pg38263@alunos.uminho.pt
- Project Link - https://github.com/deBeir/GeneAnalysis
