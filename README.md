# Machaon

This repository contains an implementation for the method presented in the  paper "Identifying and   
profiling structural similarities between Spike of SARS-CoV-2 and other viral or host proteins with  
Machaon".

<u>Please consult this time-saving manual before you use Machaon</u>. It contains an in-depth explanation    
about installing, setting up and using this method.

## System Requirements

The target system for Machaon's development is Ubuntu 20.4. Machaon has limited functionality  
on Windows and MacOS. Some post-processing modules utilize TM-Align and DSSP which are not     
cross-platform implementations. DSSP data might also be used for setting the targets of constrained    
comparisons, which is Machaon's default behaviour.    

The recommended ways to use Machaon is either by working inside a Docker container or a Singularity   
container or by working in an Ubuntu 20.4 environment with Anaconda (see instructions in the 'Installation'   
section below). On Windows, you could try WSL in order to get access to a UNIX environment (not tested):   
https://docs.microsoft.com/en-us/windows/wsl/install

Machaon is an I/O (input/output) intensive implementation and the performance is depended on the    
storage hardware and the storage optimizations of the host operating and file systems. For every   
PDB file that is analyzed, there is a corresponding set of serialized data objects in the form of   
binary files (pickle Python package) which hold the necessary data for the calculation of each    
metric. NVMe storage is highly recommended.

Machaon is a multi-core CPU application with moderate demands on RAM memory only for  
post-processing  and target setup for constrained comparisons due to the required alignments   
(especially alignments in parallel).


## Repository contents

- assess: this folder contains scripts for Machaon's benchmarking, evaluation and assessment
- config: configuration files
- docs: It contains programming-related documentation and diagrams.
   - docs/classes: Extensive API documentation for all the classes of this implementation.   
  Each class has a dedicated HTML file with thorough description.
- setup: Scripts for downloading and preparing some (optional) related data sources. 
- src: source code 
- test: It contains an integrity test with testing data and expected outputs.
- docker-compose.yml : A file used by Docker Compose tool.
- Dockerfile:  A file with the commands needed to set up Machaon in a Docker container.
- environment.yml: A file used by Anaconda Python package manager.
- LICENSE.md: The license of this implementation.
- README.md: Machaon's manual (the one you are reading). 
- Singularity: A file used to set up a Singularity container.

## Setup instructions

### Local data sources

Enrichment and meta-analysis stages rely on external data sources.  There are fallbacks in place for  
some of them (webservice calls) but it is strongly recommended utilizing the available static resources.  
This will minimize network activity, greatly speed up the process and protect the respective third party   
web services from burden. Be sure to have enough available disk space (at least 30GB) for the initial  
downloads (at least 12GB after the preparation).


<b>Note</b>: You can use the <b>'noThirdPartyData'</b> flag in the configuration, ending up only with the comparison   
results. This mode does not require the set up of local data sources or other external data access. The metrics  
<b>do not rely on external information </b> apart from the PDB file. Therefore, you only need to collect a set of   
PDB files to compare with your PDB of choice . However, you will miss enrichment and gene ID-based filtering   
of the results along with the functionality of the evaluation, meta-analysis, presentation modules.   
Also, you will not able to perform the domain scanning since it requires the residue positions of the domains   
(information found in UniProt data).

Choose a folder that will be the root data & cache folder of Machaon and <b>copy</b> there the .sh files located   
in the setup folder. You can use symbolic links if you need to have some resources in separate locations  
(https://en.wikipedia.org/wiki/Symbolic_link). Make sure the scripts have adequate execution permissions:  
```chmod 770 *.sh```

#### PDB files:

There are two ways that you can obtain multiple PDB files:

 - https://www.rcsb.org/downloads   
The downloaded files need to be decompressed and renamed (execute in the downloaded files' folder):  
```ls *.gz | parallel gunzip```   
```for f in *.ent; do mv -- "$f" "${f%.ent}.pdb"; done```   
<br/>
 - (Unix or MacOS only) https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script  
The downloaded files need to be decompressed (execute in the downloaded files' folder):   
```ls *.gz | parallel gunzip```

Of course, you can use RCSB search and retrieve relevant PDB IDs by a query of choice.

<b>Note</b>: PDB files from AlphaFold's predictions are <b> fully </b> supported. You can download them from here:  
https://alphafold.ebi.ac.uk/download

<b>Important:</b> Avoid underscores in custom PDB filenames. For example, in Ubuntu you can run:   
```rename.ul '_' '' *.pdb``` and remove an underscores from every filename in the folder.    


#### RefSeq:  

- Execute ``` ./download_refseq_resources.sh ```   
If there are any errors during the downloads, you could try to run the script a while
later (https://www.biostars.org/p/493656).
- Execute ``` ./download_refseq_resources.sh ``` again for a final verification of the
downloaded files' integrity and then execute:   
``` ./prepare_refseq_resources.sh ``` 

#### Uniprot mapping:

- It is recommended to use a dedicated FTP transferring program than a browser for the following large  
downloads (e.g. FileZilla: https://filezilla-project.org/download.php)
- Visit the following directory : https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping
- Download the following files: idmapping_selected.tab.gz, idmapping.dat.gz (Be sure to have enough space for the downloads)
- Execute ``` ./prepare_uniprot_resources.sh ``` 


### Installation (Containers)


#### Docker

If you are going to use Docker, you only need to specify your data storage in docker-compose.yml file:    
```- MY_BIG_STORAGE_PATH:/opt/storage```    
(replace MY_BIG_STORAGE_PATH with your path of choice)

and run the following command to build and launch a Machaon-ready container:    
```sudo docker-compose up -d``` 

You can enter into the container and start working:   
```sudo docker exec -it <container's name> bash```     
```conda activate machaon```   
```cd /opt/src```   
```python run.py -h```

The folder with the configurations (config) is the shared between the host system    
and  container for ease of use (you can read and edit configuration files outside of   
the container).

Alternatively, if you plan to run it in a Cloud VM instance, you need to modify the    
Docker configurations:

- docker-compose.yml: Set your mounts accordingly (or remove the volume directive)
- Dockerfile: Add the following line before WORKDIR command:   
```ADD ./config /opt/machaon/config```   
<br/>

#### Singularity

These are the instructions for creating a container with Singularity (https://sylabs.io/docs/):     

- Download the latest version from here: https://github.com/sylabs/singularity/releases
- Execute:   
```singularity build --fakeroot machaon.sif Singularity```    
```singularity run machaon.sif```   
```conda activate machaon```    
```cd /opt/src```  
```python run.py -h```
<br/><br/>  

### Manual Installation

This section is a walkthrough for manual installation (<u>please also check Dockerfile</u>, it contains all    
needed commands but it is recommended to execute them separately). 

#### Modified TM-Align compilation 

This well-established method is used for 3D similarity computation by the evaluation module.  
Machaon can run without the presence of this executable but you will miss the 3D similarity    
evaluation of the final candidates in the Machaon's results.

According to the original documentation, TM-Align is compiled as:   
```g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp```   
(You might need to install g++ first: ```sudo apt-get install build-essential``` )   
MacOS users should omit '-static' option.
For more, you can check: https://zhanggroup.org/TM-align    
     
<br/>

#### DSSP 

This well-known method is used for protein secondary structure assignment, employed in constrained   
search mode and the Gene Ontology meta-analysis process of Machaon. Alternatively, you could use   
protein or hydrophobicity-focused sequences that do not require this program otherwise Machaon  
will use STRIDE instead (see next section).

Below are the steps for the compilation of DSSP 4.0 in <b>Ubuntu 20.4</b>:

- CMake:  
```sudo apt-get install cmake```    


- Boost:  
```sudo apt-get install libboost-all-dev```  
<br/> 
   For Ubuntu versions lower than 20.04, you need to install Boost from source if your latest version is lower than 1.70: <br/><br/>  
   - Remove previous Boost version:      
```apt remove 'libboost.*-dev'```  
   - Download and extract the latest version from: https://www.boost.org/  (greater than 1.70)  
   - Install:  
```chmod +x bootstrap.sh```  
```./boostrap.sh```  
```sudo ./b2 link=static install```  


- BZIP2:  
```sudo apt-get install libbz2-dev```


- cifpp:
Make sure you have cmake (```sudo apt install cmake ```) and follow the instructions in https://github.com/PDB-REDO/libcifpp  
You might need also to install this before: https://github.com/mhekkel/mrc (https://github.com/PDB-REDO/dssp/issues/4)  
For Ubuntu 18.04 you also need to install these first of all:  
```sudo add-apt-repository ppa:ubuntu-toolchain-r/test```  
```sudo apt update```  
```sudo apt install gcc-9 g++-9```  
```export CC=/usr/bin/gcc-9```  
```export CXX=/usr/bin/g++-9```  


- DSSP: Please follow the instructions in  https://github.com/PDB-REDO/dssp

<b>Note:</b>There are also other options to obtain DSSP files without setting up the program: https://swift.cmbi.umcn.nl/gv/dssp/  
In that case, you should add them in a folder named 'dssp_cache' located in your specified root data & cache folder   
('rootDisk' parameter, more in 'Execution' section) . 
<br/><br/>

#### STRIDE

STRIDE is an established method for determining the protein secondary structure from PDB files.
It is used as a fallback solution for custom PDB files that do not fully follow the standard PDB 
format and lack annotations. Please follow the instructions in  http://webclu.bio.wzw.tum.de/stride/

<br/>
After the compilations, you have to copy the mkdssp, stride, TM-Align executables  
into the directory of Machaon and give them the required execute permissions:   

```cd machaon/src```  
```cp /home/<user>/.local/bin/mkdssp .```  
```cp /home/<user>/.local/bin/stride .```  
```chmod 770 mkdssp ```  
```chmod 770 stride ```  
```chmod 770 TMalign ```  
<br/><br/>
 

#### Required system libraries:

You need the poppler library in order to export the figures in the EPS format
with Python plotly library:  
```sudo apt-get install  libpoppler-cpp-dev```
This a graphics related library for Open3D:
```sudo apt-get install  libgl1-mesa-dev```
<br/><br/>

#### Python environment:

An environment setup of Anaconda Python distribution is needed : https://www.anaconda.com

This distribution allows easy setup of all the requisites for Machaon.

Once you have an operational Anaconda-enabled terminal, move into the setup folder and execute  
the following command to install all the required packages:  
```conda env create -f environment.yml```  
<br/>

### Testing your installation:

Run the test script in the /test folder:
``` python integrity_test.py ```   
If there are no differences reported at the end, than your installation should be successful.

## Execution

At first, you need to activate the previously installed environment in an Anaconda-enabled terminal:   
```conda activate machaon```

### Quick start:  

Execute the following script which is located in the src folder:  ``` run.py -h```  
This will display all the available options and their descriptions.


### Batch jobs (recommended):
Edit <b>config.yaml</b> file in the src folder and run <b> batch_run.py</b>. Below is an example entry with the default  
values. You could copy it and modify it according to your needs. Configurations with "ignore : True" field  
are ignored. You could also consult with the example configurations used for the Spike protein of SARS-CoV-2.
<pre>
  - rootDisk: "" 
    referencePDBID: ""
    overridePDBID: ""
    referenceChainID: ""
    referenceGeneID: ""
    referenceSequenceLength: 0
    comparisonMode: ""
    pdbDatasetPath: ""
    outputPath: ""
    excludedOrganisms: []
    excludedGeneNames: []
    excludedPDBIDs: []
    isReferenceViral: False
    GOProperty: ""
    GOTargetProperties: []
    GOSearch: ""
    GOAlignmentLevel: "secondary"
    noThirdPartyData: False
    pdbValidation: False
    GOAnalysisOnly: False 
</pre>
All the options are presented below:
<pre>
'rootDisk': This will also be the caching location for the extracted features.
'referencePDBID': Choose the reference PDB IDs (1 search per reference)
'overridePDBID': Override the reference PDBID for Uniprot ID retrieval (for renamed reference PDB files, e.g. 6VXX_processed.pdb)
'referenceChainID': Choose the chain of the reference PDB
'referenceGeneID': Provide the gene id (Entrez) of the reference PDB
'referenceSequenceLength':  Provide the protein sequence length of the reference protein
'comparisonMode':  Choose 'whole', 'domain' or 'segment'
'alignmentLevel': Choose 'primary', 'secondary', 'mixed', 'hydrophobicity'. Default is 'mixed'. (Only from segment scans)
'pdbDatasetPath': Relative path for PDB data folder
'outputPath': The location of the outputs (can be relative or full path)
'excludedOrganisms': Filtering out structures originating from the same organism as the reference one
'excludedGeneNames': Filtering out structures originating from the same gene as the reference one
'excludedPDBIDs': Exclude PDB IDs
'isReferenceViral': Meta-analysis skips the search in viral genome data for the reference, if it is not a viral protein
'GOProperty': Choose a property type for analysis: 'biologicalProcess', 'molecularFunction', 'cellularComponent'
'GOTargetProperties': Choose properties for analysis
'GOSearch': Choose a term to be searched in all available GO Terms belonging to the results e.g. 'ubiquit' (could be a stem of a word)
'GOAlignmentLevel': Choose target alignment level : ['primary', 'secondary', 'mixed', 'hydrophobicity'. Default is 'secondary']
'noThirdPartyData': Do not use external local or online resources. PDB data only.
'GOAnalysisOnly': Perform only GO Meta-analysis (for completed searches).
'pdbValidation': Validation for PDB files. Every file assessed as invalid is skipped from the search (very strict and slow). 
'ignore': If set to True, the configuration will be ignored. Useful for storing previous job details.
'*whatever*': You can include fields of your own, like tags or notes (e.g. 'date' : 14-4-2003). These are not considered by the program. 
</pre>

#### Constrained mode (segments):
Constrained search on segments requires also preset about the reference segment. This is set in  
<b>segments.yaml</b> file in the src folder. Below is an empty template entry.  You could also consult with  
the example segment definitions used for the Spike protein of SARS-CoV-2.

<pre>
- referencePDBChain: ""
  residues: []
  residueRanges: ""
  known: False
</pre>

All the options are presented below:
<pre>
'referencePDBChain': The reference PDB ID and chain ID separated by a dot, &lt;PDB ID&gt;.&lt;CHAIN ID&gt; e.g. ""6VXX.A"
'residues': List of residue positions (one-based indexing), e.g. [1, 2, 3, 4, 5]
'residueRanges': Range definitions separated by comma, e.g. '1-50,70-78'
'known': Select True if the segment belongs to a known site like a binding site (considered by GO Meta-analysis module).
'ignore': If set to True, the configuration will be ignored. Useful for storing past segment presets.
'*whatever*': You can include fields of your own, like tags or notes (e.g. 'doi' : '123.3456/1234.123'). These are not considered by the program. 

Note: 'residues' and 'residueRanges' definitions are combined, e.g. [12, 15, 59]  
and '13-40, 47-52' would result to the selection of residue positions from 12 to 40, 
from 47 to 52 and 59 (duplicate definitions are removed).
</pre>

## Data directory structures

#### Output folder structure
<pre>
    (a user-specified output folder)
    |
    |__ (a user-specified top directory name)
        |
        |__metrics/ (directory for the computed metrics for all structures in the dataset)
        |
        |__candidates/ (directory for the selected final set of candidate entries, 
        |               the final report is saved here [HTML file])
        |
        |__plots/ (directory for plots regarding the final set)
        |
        |__go/ (directory for GO meta-analysis and related visualizations)
</pre>

<b>Note for constrained mode search on segments</b>:The corresponding output files contain a suffix   
"site&lt;segment index&gt;" that signify the results for a particular segment. The index comes from the   
configuration order. In the "metrics" folder, there is a "*_site&lt;segment index&gt;-parts.csv" file that contains  
the contiguous parts of the segment as determined by the method.


#### Root folder (source data & cache), full structure

<pre>
    (a user-specified <b>root</b> folder)
    |
    |--DATA_&lt;PDB directory name&gt;_&lt;whole or domain&gt;/ 
    |  (directory for storing the extracted features of a PDB directory)
    |
    |--domains/ (directory for caching domain information by UniProt online requests)
    |
    |--dssp_cache/ (directory for caching DSSP results)
    |
    |--enrichment/ (directory for caching data enrichment of PDB chain entries)
    |
    |__entrez/ (cache directory for NCBI Entrez online requests)
    |
    |--pdbinfo/ (directory for caching extracted PDB meta-data)
    |
    |--prot_sec/ (directory for caching PDB sequence/secondary structure data)
    |
    |__refseq/ (RefSeq resources directory)
    |
    |--rcsbenrich/ (cache directory for RCSB enrichment data) 
    |
    |--(user created PDB folders, <b>each folder corresponds to a target dataset for a search</b>)
    |
    |__idmapping_selected.tab.gz (UniProt idmapping resources)
</pre>

There is also a cache file that is generated besides the scripts in src folder (go_cache.csv) that holds  
Gene Ontology data.

## Output format

The outputs are human interpretable CSV files with headers:

- metrics directory has comma separated CSV files
- candidates directory has tab separated CSV files
- outputs of constrained searches include columns with serialized list contents which can be parsed with eval()

## Special Cases

- If you want to compare a polymer as a whole structure you could use pdb-tools :
https://github.com/haddocking/pdb-tools   
and combine multiple chains to one. You should remove any pre-computed features of the old PDB  
(*_angles.pkl, *_distances.pkl, *_triangles.pkl) and the original PDB from the dataset (you could  
keep these files in a separate location as back up). You need to decide which original &lt;PDB ID&gt; and  
&lt;PDB chain ID&gt; you will use as a reference for the third-party resources.

- In case you encounter warnings about empty chain identifiers or missing chains, use pdb_chain  
command from pdb-tools: ```pdb_chain -A no_chains.pdb > corrected.pdb``` to put a dummy identifier
to a problematic PDB file.

- Secondary structure data cannot be extracted from PDBs that lack experimental information so you may have to 
change the target alignment level to primary or hydrophobicity (recommended) for constrained mode search on 
segments (default is 'mixed') or GO metanalysis (default is 2D). 





### Trivia  
  
https://en.wikipedia.org/wiki/Machaon_(mythology)
