import configparser  
import json 
import os
import sys
from Bio import Entrez
import pandas as pd
from tqdm import tqdm
from pathlib import Path


def read_config():
    config = configparser.ConfigParser()
    config.read('config.ini')
    # Access variables directly from any section
    email = config.get('Pubmed', 'email')
    api_key = config.get('Pubmed', 'api_key')
    retmax = int(config.get('Pubmed', 'retmax'))    
    output_file = config.get('Pubmed', 'main_output')
    term_list = json.loads(config.get('Pubmed', 'term_list'))  # This reads the actual list
    output_folder = config.get('Pubmed', 'output_folder')
    list_output_folder = config.get('Pubmed', 'list_output_folder') 
    PM_output_folder = os.path.join(output_folder, list_output_folder)
    os.makedirs(PM_output_folder, exist_ok=True)
    return  email, api_key, retmax, output_file, term_list, output_folder, PM_output_folder

# Function to fetch article details and extract DOI
def fetch_article_data(pmid):
    handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
    rec = Entrez.read(handle)
    handle.close()

    article = rec["PubmedArticle"][0]["MedlineCitation"]["Article"]

    # Extract DOI from ArticleIdList
    doi = None
    try:
        article_ids = rec["PubmedArticle"][0]["PubmedData"]["ArticleIdList"]
        for aid in article_ids:
            if aid.attributes["IdType"].lower() == "doi":
                doi = str(aid)
                break
    except Exception:
        doi = None

    # Format DOI for Excel
    doi_excel = f'=HYPERLINK("https://doi.org/{doi}", "{doi}")' if doi else ""

    # Extract authors
    authors = []
    for author in article.get("AuthorList", []):
        last = author.get("LastName", "")
        first = author.get("ForeName", "")
        if last and first:
            authors.append(f"{last} {first}")
        elif last:
            authors.append(last)
    authors_str = ", ".join(authors)

    # Extract other fields
    pmid = rec["PubmedArticle"][0]["MedlineCitation"].get("PMID", "")
    title = article.get("ArticleTitle", "")
    journal = article.get("Journal", {}).get("Title", "")
    issn = article.get("Journal", {}).get("ISSN", "")
    year = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {}).get("Year", "")

    # Abstract extraction: could be a list or string
    abstract = ""
    if "Abstract" in article:
        abs_obj = article["Abstract"].get("AbstractText", "")
        if isinstance(abs_obj, list):
            abstract = " ".join([str(a) for a in abs_obj])
        else:
            abstract = str(abs_obj)

    return {
        "PMID": pmid,
        "Title": title,
        "Authors": authors_str,
        "Journal": journal,
        "ISSN": issn,
        "Year": year,
        "DOI": doi_excel,
        "Abstract": abstract
    }
   
#Merge all CSVs and drop duplicate PMIDs
def merge_csv_folder(PM_output_folder: str, output_name: str = "01_OUTPUT.csv", output_folder: str = "./Outputs/") -> Path:
    folder = Path(PM_output_folder)
    if not folder.is_dir():
        raise ValueError(f"Not a directory: {folder}")
    csv_files = sorted(folder.glob("*_PMoutput.csv"))
    if not csv_files:
        raise ValueError(f"No CSV files found in {folder}")

    dfs = []
    for p in csv_files:
        try:
            df = pd.read_csv(p, dtype=str, low_memory=False, on_bad_lines="skip")
            # normalize column names (strip)
            df.columns = [c.strip() if isinstance(c, str) else c for c in df.columns]
            dfs.append(df)
        except Exception as e:
            print(f"[WARN] Failed to read {p.name}: {e}", file= sys.stderr)

    if not dfs:
        raise RuntimeError("No readable CSVs found")

    # align columns by union (pandas.concat will fill NaN for missing columns)
    merged = pd.concat(dfs, ignore_index=True, sort=False)

    # detect PMID column
    before = len(merged)
    merged = merged.drop_duplicates(subset=['PMID'], keep="first")
    after = len(merged)
    print(f"Dropped {before - after} duplicate rows by column 'PMID'")

    out_path = Path(output_folder) / output_name
    merged.to_csv(out_path, index=False)
    print(f"Saved merged table to: {out_path}")
    return out_path


def Literature_search(term_list, email, api_key, output_folder, output_file, PM_output_folder, retmax):
    print("Configuration loaded successfully.")
    # Set your email (required by NCBI)
    Entrez.email = email
    Entrez.api_key = api_key
    n = 0
    #term_list = term_list[:5]
    # Loop over each tissue and run the search
    for tissue in term_list:
        n += 1
        search_term = f"Phosphoproteom* AND {tissue} AND (mass AND spectrometry)"
        print(f"Searching for: {search_term}, number {n} out of {len(term_list)}")
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]
        # Step 2: Fetch and parse data with progress bar
        data = []
        for pmid in tqdm(id_list, desc=f"Fetching articles for {tissue}"):
            try:
                article_data = fetch_article_data(pmid)
                data.append(article_data)
            except Exception as e:
                print(f"Error fetching PMID {pmid}: {e}")

        tissue_output_file = f"{tissue.replace(' ', '_')}_PMoutput.csv"
        df = pd.DataFrame(data)
        output_path = os.path.join(PM_output_folder, tissue_output_file)
        if df.empty:
            print(f"No data found for {tissue}")
        else:
            df.to_csv(output_path, index=False)
            print(f"Saved {len(df.index)} results to {tissue_output_file}")
            #_______________________________________________________________________
    #H_index match

    sjr_df = pd.read_csv("Inputs/scimagojr 2024.csv", sep=';')
    # Prepare SJR ISSNs and H index
    sjr_df = sjr_df.dropna(subset=['Issn'])
    sjr_exploded = sjr_df.assign(ISSN_split=sjr_df['Issn'].str.split(',')).explode('ISSN_split')
    sjr_exploded['ISSN_split'] = sjr_exploded['ISSN_split'].str.replace('-', '', regex=False).str.upper().str.strip()
    sjr_hindex = sjr_exploded[['ISSN_split', 'H index']].drop_duplicates()

    # Find all CSV files ending with _PMoutput.csv in the output folder
    import glob
    csv_files = glob.glob(os.path.join(PM_output_folder, "*_PMoutput.csv"))
    
    for csv_file_path in csv_files:
        tissue_output_file = os.path.basename(csv_file_path)
        output_path = csv_file_path

        # Load data
        try:
            pubmed_df = pd.read_csv(output_path)
            if pubmed_df.empty:
                print(f"No data found in {tissue_output_file}")
                continue
            
            # Clean and normalize ISSNs in pubmed_df
            if 'ISSN' in pubmed_df.columns:
                pubmed_df['ISSN_clean'] = pubmed_df['ISSN'].str.extract(r'(\d{4}-\d{3}[0-9Xx])')
                pubmed_df['ISSN_clean'] = pubmed_df['ISSN_clean'].str.replace('-', '', regex=False).str.upper().str.strip()
            else:
                print(f"Warning: No 'ISSN' column found in {tissue_output_file}")
                pubmed_df['ISSN_clean'] = None

            # Merge H index into pubmed_df; use suffixes to avoid immediate collision
            merged_df = pd.merge(pubmed_df, sjr_hindex, left_on='ISSN_clean', right_on='ISSN_split', how='left', suffixes=('_left','_right'))

            # If merge produced paired columns like 'H index_left' and 'H index_right', coalesce them into 'H index'
            cols_to_drop = []
            for col in list(merged_df.columns):
                if col.endswith('_left'):
                    base = col[:-5]
                    right_col = base + '_right'
                    if right_col in merged_df.columns:
                        merged_df[base] = merged_df[col].combine_first(merged_df[right_col])
                        cols_to_drop.extend([col, right_col])
            if cols_to_drop:
                merged_df.drop(columns=cols_to_drop, inplace=True)

            # Drop the extra merge column
            merged_df.drop(columns=['ISSN_split'], inplace=True)
            merged_df.drop(columns=['Unnamed: 0'], inplace=True, errors='ignore')
            merged_df.set_index('PMID', inplace=True)
            merged_df.to_csv(f"{output_path}")
            print(f"Processed {tissue_output_file} successfully.")

        except pd.errors.EmptyDataError:
            print(f"Empty data error for {tissue_output_file}")
            continue
    merge_csv_folder(PM_output_folder=PM_output_folder, output_name=output_file, output_folder=output_folder)
            

#if __name__ == "__main__":
#    email, api_key, retmax, output_file, term_list, output_folder, PM_output_folder = read_config()
#    Literature_search(term_list, email, api_key, output_folder, output_file, PM_output_folder, retmax)
      
