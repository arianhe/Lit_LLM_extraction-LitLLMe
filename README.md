# Lit_LLM_extraction-LitLLMe
Tool for extensive literature research across disease types. Uses kinase inhibitor lists to automate PubMed searches, maps results to ProteomeXchange and PRIDE, and extracts metadata and PDFs. These serve as an input for further LLM extraction.

## Before execution of the pipeline

Setup of environment
```
conda env create -f environment.yml
conda activate jupyter_env
```
Download the appropriate CUDA compatible llama-cpp from GitHub. 

For Linux inside the console:  
```
 CMAKE_ARGS="-DGGML_CUDA=on -DLLAMA_CUDA=on" pip install git+https://github.com/abetlen/llama-cpp-python.git --no-cache-dir
 ```

Execute the LLM download.


Additionally you will need to specify an email adress and your NCBI API-key in the config file in order to speed up the process.
 
Start the literature search by selecting your preferred search term in `Supp_Literatur_table_OpenTargetInt.py` and/or use the list in config.ini to cycle through a list of terms.  
Review for query grammar: https://biopython.org/docs/1.84/Tutorial/chapter_entrez.html and https://www.ncbi.nlm.nih.gov/books/NBK25499/  

## Pipeline  
 * 00_Installing_LLM.py
	 * Installs the model capybarahermes-2.5-mistral-7b.Q5_K_M in a seperate model folder. As the file is .gguf an installation via this script is recommended (browser download can damage the file due to compression)
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 
 * Supp_Literatur_table_OpenTargetInt.py
   * *no need to execute, contains functions for pubmed query*
 * A_Download_OpenTarget.py
   * Cycles through terms from `config.ini` (`term_list`), queries OpenTarget one at a time
   * Selects the first hit and extracts **unique kinase-related drugs**
   * For each drug, runs a PubMed search using:
    
    ```python
    search_term = f"Phosphoproteom* AND {tissue} AND (mass AND spectrometry)"
    ```
    
   * Outputs:
     - Individual `_PMoutput.csv` files with: **PMID, Title, Authors, Journal, ISSN, Year, DOI, Abstract**
     - Merged file: `01_OUTPUT.csv`
     - Folder: `Outputs/Pubmed_outputs/` for single results
   * Maps papers to **journal h-index** using the merged CSV.
 * B_1_ProtXchange.py  
     * This script will download all of the metadata regarding every project on PRIDE as a json (filesize ~ 240 MB), fix parsing errors. Then it will extract 'accession', 'title', 'doi', 'organisms', 'tissues', 'diseases', 'instruments', 'projectDescription', 'sampleProcessingProtocol', 'experimentType', 'quantificationMethod' from "all_pride_projects_fixed.json".  This is modifiable in 02_PRIDE.py.
	 * The resulting output is saved into "all_pride_projects_extracted.csv".  
     * This information is still provided by PRIDE and the extraction is modifiable. Check the structure of "all_pride_projects_fixed.json".   
     * Afterwards it will map the results of the literature search "01_OUTPUTS.csv" to extracted PRIDE projects by PMID (Pubmed id) and merge into one dataframe.
     * As a fallback the metadata from ProteomeXchange is downloaded too as `proteomexchange_all` and mapped to PMIDs that could not be mapped to PRIDE.
     * The merged output is accessible as "02_OUTPUT.csv" while single matches are accessible as "*_pmid_matches.csv" in Outputs/Pubmed_outputs.
     * Configuration whether only PRIDE matches shall be kept, can be made in the config file at `drop_non_PRIDE`.
 * B_2_PDF_download.py
 	 * This script will read any PMID in the PMID column of a file specified in config.ini and try to download it as a pdf (PMID.pdf) to Outputs/Papers.
     * Then it will read the filenames and check for any missing PMIDs and save them to missing_pmids.txt, these require often institutional authentication.  
 * C_LLM_Extraction.py 
     * This is the heart of the extraction pipeline. This script creates a new column (later removed) consisting of the paper abstract and the sampleProcessingProtocol (seperated by ' ') and uses this as an input to extract cell lines, targeted kinases and the used drug.   
     * As an additional input to the prompt it uses a list of cell lines ("cancer_cell_lines_simplified.csv"), kinases and their inhibitors retrieved from A_Download_OpenTarget.py or ("Kinase_Inhibitors.csv", sep =";" as a *fallback*).  
     * **IMPORTANT:** The environment the script is executed in _should_ have **access** to the **GPU** due to runtime issues. 
     * The output of a single input file is saved as "03.OUTPUT.csv" whereas the output of the PDF extraction is saved as "*_extracted.csv" in the `Papers` directory.  

Have fun and happy coding! :hammer_and_wrench:
