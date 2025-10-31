import subprocess
import json
import sys
import pandas as pd
import re
import os
import configparser
from pathlib import Path
import requests
import xml.etree.ElementTree as ET
import ijson
import csv
from tqdm import tqdm


def read_config():
    config = configparser.ConfigParser()
    config.read('config.ini')
    return config

config = read_config()
output_file = config.get('Pubmed', 'main_output')
output_folder = config.get('Pubmed', 'output_folder')
list_output_folder = config.get('Pubmed', 'list_output_folder')
#PRIDE
list_of_inputs = config.getboolean('Input_paths', 'list_of_inputs')
folder = config.get('PRIDE', 'folder')
metadata_file = config.get('PRIDE', 'metadata_file')
fixed_json_file = config.get('PRIDE', 'fixed_json_file')
output_csv_file = config.get('PRIDE', 'output_csv_file')
Update_PRIDE = config.getboolean('PRIDE', 'Update_PRIDE')
drop_non_PRIDE = config.getboolean('PRIDE', 'drop_non_PRIDE')
main_output = config.get('PRIDE', 'main_output')

PM_output_folder = os.path.join(output_folder, list_output_folder)
os.makedirs(folder, exist_ok=True)
metadata_file = os.path.join(folder, metadata_file)
fixed_json_file = os.path.join(folder, fixed_json_file)
output_csv_file = os.path.join(folder, output_csv_file)

def fix_strings_control_chars(text: str) -> str:
    allowed = {9, 10, 13}
    def clean_char(c):
        if ord(c) < 32 and ord(c) not in allowed:
            return ' '
        return c
    return ''.join(clean_char(c) for c in text)

def fix_unescaped_newlines_in_json_strings(text: str) -> str:
    string_regex = re.compile(r'"(.*?)(?<!\\)"', re.DOTALL)
    def replacer(match):
        s = match.group(0)
        inner = s[1:-1]
        inner_fixed = inner.replace('\n', '\\n').replace('\r', '\\r')
        return f'"{inner_fixed}"'
    return string_regex.sub(replacer, text)

def Download():
    if Update_PRIDE == False:
        if os.path.exists(metadata_file) or os.path.exists(fixed_json_file):
            print(f"{metadata_file} or {fixed_json_file} already exists. Skipping execution.")
            return
    else:
        print("Update_PRIDE is set to True. Proceeding with download and processing.")

    # Step 1: Download all PRIDE project metadata via pridepy CLI
    print("Downloading all PRIDE metadata...")
    subprocess.run(["pridepy", "stream-projects-metadata", "-o", metadata_file], check=True)

    print("Reading raw JSON...")
    with open(metadata_file, "r", encoding="utf-8", errors='replace') as f:
        raw_text = f.read()

    print("Fixing control characters...")
    fixed_text = fix_strings_control_chars(raw_text)

    print("Fixing unescaped newlines in JSON strings...")
    fixed_text = fix_unescaped_newlines_in_json_strings(fixed_text)

    print(f"Writing cleaned JSON to {fixed_json_file} ...")
    with open(fixed_json_file, "w", encoding="utf-8") as f:
        f.write(fixed_text)

    print("Done: Cleaned JSON written to", fixed_json_file)

def download_all_proteomexchange_metadata(output_file: str = "proteomexchange_all.xml"):
    """Download metadata from ProteomeXchange central registry (all repos)."""
    url = "http://proteomecentral.proteomexchange.org/cgi/GetDataset"
    params = {"outputMode": "XML", "test": "no"}
    
    print("[ProteomeXchange] Downloading all repository metadata...")
    try:
        resp = requests.get(url, params=params, timeout=120)
        resp.raise_for_status()
        Path(output_file).write_bytes(resp.content)
        print(f"[ProteomeXchange] Metadata saved to {output_file}")
        return True
    except Exception as e:
        print(f"[ProteomeXchange][ERROR] Download failed: {e}")
        return False

def parse_proteomexchange_xml(xml_file: str):
    """Parse ProteomeXchange XML and extract project records with PMIDs."""
    tree = ET.parse(xml_file)
    root = tree.getroot()
    projects = []
    
    for dataset in root.findall(".//Dataset_Identifier"):
        dataset_id = dataset.get("id", "")
        
        accession = ""
        if dataset_id:
            pxd_match = re.search(r'(PXD\d+)', dataset_id)
            if pxd_match:
                accession = pxd_match.group(1)
            else:
                accession = re.sub(r'<[^>]+>', '', dataset_id).strip()
        
        title_elem = dataset.find("Title")
        title = title_elem.text.strip() if title_elem is not None and title_elem.text is not None else ""
        
        repos_elem = dataset.find("Repos")
        repos = repos_elem.text.strip() if repos_elem is not None and repos_elem.text is not None else ""
        
        pmid = extract_pmid_from_proteomexchange_publication(dataset)
        
        projects.append({
            "accession": accession,
            "title": title,
            "pmid": pmid,
            "repos": repos
        })
    
    return projects

def extract_pmid_from_proteomexchange_publication(dataset):
    """Extract PMID from ProteomeXchange Publication field."""
    pub_elem = dataset.find("Publication")
    if pub_elem is None:
        return ""
    
    pub_text = ""
    if pub_elem.text:
        pub_text += pub_elem.text
    
    for child in pub_elem:
        if child.get('href'):
            pub_text += " " + child.get('href')
        if child.text:
            pub_text += " " + child.text
        if child.tail:
            pub_text += " " + child.tail
    
    pub_text = pub_text.strip()
    if not pub_text:
        return ""
    
    try:
        pmid_pattern = r'https://www\.ncbi\.nlm\.nih\.gov/pubmed/(\d+)'
        match = re.search(pmid_pattern, pub_text)
        if match:
            return match.group(1)
        
        pmid_alt_pattern = r'PMID[:\s]*(\d+)'
        alt_match = re.search(pmid_alt_pattern, pub_text, re.IGNORECASE)
        if alt_match:
            return alt_match.group(1)
        
        pmid_numeric_pattern = r'\b(\d{6,})\b'
        numeric_matches = re.findall(pmid_numeric_pattern, pub_text)
        if numeric_matches:
            return numeric_matches[0]
            
    except Exception:
        return ""
    
    return ""

def create_proteomexchange_dataframe(xml_file: str = "proteomexchange_all.xml"):
    """Create a pandas DataFrame from ProteomeXchange XML data."""
    try:
        projects = parse_proteomexchange_xml(xml_file)
        
        if not projects:
            return pd.DataFrame()
        
        df = pd.DataFrame(projects)
        df['PMID_clean'] = df['pmid'].apply(clean_pmid)
        
        df = df[df['PMID_clean'].notna()]
        df = df.drop_duplicates(subset=['PMID_clean'])
        
        print(f"ProteomeXchange: {len(df)} records with valid PMIDs")
        return df
        
    except Exception as e:
        print(f"Error processing ProteomeXchange data: {e}")
        return pd.DataFrame()

def extract_pmid_from_references(refs):
    if isinstance(refs, list):
        for ref in refs:
            if 'pubmedID' in ref and ref['pubmedID']:
                return str(ref['pubmedID']).strip()
    return ''

def extract_cv_values(entries):
    """Extracts 'name' from a list of CV term dicts."""
    if isinstance(entries, list):
        return "; ".join(entry.get("name", "") for entry in entries if entry.get("name"))
    return ''

def process_projects(input_json_path, output_csv_path):
    with open(input_json_path, 'r', encoding='utf-8') as f:
        projects_iter = ijson.items(f, 'item')

        with open(output_csv_path, 'w', encoding='utf-8', newline='') as csvfile:
            fieldnames = [
                'accession', 'title', 'pmid',
                'organisms', 'tissues', 'diseases',
                'instruments', 'projectDescription',
                'sampleProcessingProtocol', 'experimentType',
                'quantificationMethod'
            ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            for project in tqdm(projects_iter, desc="Processing projects"):
                try:
                    accession = project.get('accession', '')
                    title = project.get('title', '')

                    pmid = extract_pmid_from_references(project.get('references', []))
                    organisms = extract_cv_values(project.get('organisms', []))
                    tissues = extract_cv_values(project.get('organismParts', []))
                    diseases = extract_cv_values(project.get('diseases', []))
                    instruments = extract_cv_values(project.get('instruments', []))
                    projectDescription = project.get('projectDescription', '')
                    sampleProcessingProtocol = project.get('sampleProcessingProtocol', '')

                    experimentTypes = project.get('experimentTypes', [])
                    experimentType = ''
                    if isinstance(experimentTypes, list) and len(experimentTypes) > 0:
                        experimentType = experimentTypes[0].get('name', '')

                    quantificationMethod = extract_cv_values(project.get('quantificationMethods', []))

                    writer.writerow({
                        'accession': accession,
                        'title': title,
                        'pmid': pmid,
                        'organisms': organisms,
                        'tissues': tissues,
                        'diseases': diseases,
                        'instruments': instruments,
                        'projectDescription': projectDescription,
                        'sampleProcessingProtocol': sampleProcessingProtocol,
                        'experimentType': experimentType,
                        'quantificationMethod': quantificationMethod
                    })

                except Exception as e:
                    print(f"⚠️ Error parsing project {accession}: {e}")

def Extract():
    if not os.path.exists(fixed_json_file):
        print(f"{fixed_json_file} does not exist. Please run the download step first.")
        return
    if os.path.exists(output_csv_file):
        print(f"{output_csv_file} already exists. Skipping extraction.")
        return
    print("Parsing cleaned JSON and writing CSV...")
    process_projects(fixed_json_file, output_csv_file)
    print("Done: Extracted metadata written to", output_csv_file)

def clean_pmid(raw_pmid):
    """Clean and standardize PMID values with validation"""
    if pd.isna(raw_pmid):
        return None
    pmid_str = str(raw_pmid).strip()
    if pmid_str == '' or pmid_str.lower() == 'nan':
        return None
    
    try:
        pmid_int = int(float(pmid_str))
        return str(pmid_int)
    except (ValueError, TypeError):
        return None

def load_and_clean_pride_data(pride_path):
    """Load and clean PRIDE data once to reuse for all tissue mappings"""
    print("Loading and cleaning PRIDE data...")
    
    if not os.path.exists(pride_path):
        print(f"⚠️ PRIDE data file not found: {pride_path}")
        print("Please run Download() and Extract() first.")
        return pd.DataFrame()
    
    pride_df = pd.read_csv(pride_path)
    pride_df['PMID_clean'] = pride_df['pmid'].apply(clean_pmid)
    pride_df = pride_df[pride_df['PMID_clean'].notna()]
    pride_df = pride_df.drop_duplicates(subset=['PMID_clean'])
    
    print(f"PRIDE data processed: {len(pride_df)} records with valid PMIDs")
    return pride_df

def map_pmids_with_cleaned_pride(pubmed_path, pride_df_cleaned, output_path):
    """Map PubMed data with pre-cleaned PRIDE data"""
    pubmed_df = pd.read_csv(pubmed_path)
   
    pubmed_df['PMID_clean'] = pubmed_df['PMID'].apply(clean_pmid)
    pubmed_df = pubmed_df[pubmed_df['PMID_clean'].notna()]
    pubmed_df = pubmed_df.drop_duplicates(subset=['PMID_clean'])
    
    if drop_non_PRIDE:
        merged_df = pd.merge(pubmed_df, pride_df_cleaned, how='inner', on='PMID_clean')
    else:
        merged_df = pd.merge(pubmed_df, pride_df_cleaned, how='left', on='PMID_clean')
    
    merged_df = merged_df[['PMID', 'Title', 'Authors', 'Journal', 'Year', 'H index', 
       'Abstract', 'ISSN_clean', 'DOI', 'PMID_clean', 'accession', 'organisms', 'tissues', 'diseases', 'instruments',
       'projectDescription', 'sampleProcessingProtocol', 'experimentType', 'quantificationMethod']]
    

    merged_df.to_csv(output_path, index=False)
    print(f"{len(merged_df)} matches of {len(pubmed_df)} search results written to {output_path}")

def map_pmids_with_proteomexchange_fallback(pubmed_path, pride_df_cleaned, proteomexchange_df, output_path):
    """Map PubMed data with PRIDE first, then ProteomeXchange for unmapped PMIDs."""
    pubmed_df = pd.read_csv(pubmed_path)
    pubmed_df['PMID_clean'] = pubmed_df['PMID'].apply(clean_pmid)
    pubmed_df = pubmed_df[pubmed_df['PMID_clean'].notna()]
    pubmed_df = pubmed_df.drop_duplicates(subset=['PMID_clean'])
    
    pride_merged = pd.merge(pubmed_df, pride_df_cleaned, how='left', on='PMID_clean')
    pride_mapped = pride_merged[pride_merged['accession'].notna()].copy()
    pride_unmapped = pride_merged[pride_merged['accession'].isna()].copy()
    
    print(f"PRIDE mapping: {len(pride_mapped)} mapped, {len(pride_unmapped)} unmapped")
    
    px_mapped = pd.DataFrame()
    px_unmapped = pride_unmapped.copy()
    
    if len(pride_unmapped) > 0 and proteomexchange_df is not None and len(proteomexchange_df) > 0:
        pride_unmapped_clean = pride_unmapped[['PMID', 'Title', 'Authors', 'Journal', 'Year', 'H index', 
                                                'Abstract', 'ISSN_clean', 'DOI', 'PMID_clean']].copy()
        
        px_merged = pd.merge(pride_unmapped_clean, proteomexchange_df, how='left', on='PMID_clean')
        px_mapped = px_merged[px_merged['accession'].notna()].copy()
        px_unmapped = px_merged[px_merged['accession'].isna()].copy()
        
        print(f"ProteomeXchange mapping: {len(px_mapped)} additional mapped")
    
    dfs_to_concat = []
    if len(pride_mapped) > 0:
        pride_mapped['mapping_source'] = 'PRIDE'
        pride_mapped['repos'] = 'PRIDE'
        dfs_to_concat.append(pride_mapped)
    if len(px_mapped) > 0:
        px_mapped['mapping_source'] = 'ProteomeXchange'
        dfs_to_concat.append(px_mapped)
    if not drop_non_PRIDE and len(px_unmapped) > 0:
        px_unmapped['mapping_source'] = 'Unmapped'
        px_unmapped['repos'] = ''
        dfs_to_concat.append(px_unmapped)
    
    if dfs_to_concat:
        final_df = pd.concat(dfs_to_concat, ignore_index=True)
        
        standard_columns = ['PMID', 'Title', 'Authors', 'Journal', 'Year', 'H index', 
                           'Abstract', 'ISSN_clean', 'DOI', 'PMID_clean', 'accession', 'mapping_source']
        
        for col in standard_columns:
            if col not in final_df.columns:
                final_df[col] = ''
        
        metadata_columns = []
        if 'organisms' in final_df.columns:
            metadata_columns.extend(['organisms', 'tissues', 'diseases', 'instruments',
                                   'projectDescription', 'sampleProcessingProtocol', 
                                   'experimentType', 'quantificationMethod'])
        if 'repos' in final_df.columns:
            metadata_columns.extend(['repos', 'title'])
        
        final_columns = standard_columns + [col for col in metadata_columns if col in final_df.columns]
        final_df = final_df[final_columns]
        
        final_df.to_csv(output_path, index=False)
        mapped_count = len(final_df[final_df['mapping_source'] != 'Unmapped'])
        print(f"{mapped_count} matches of {len(pubmed_df)} mapped to {output_path}")
        print(f"{len(final_df)} total records written to {output_path}")
        
        mapping_summary = final_df['mapping_source'].value_counts()
        print("Mapping summary:")
        for source, count in mapping_summary.items():
            print(f"  {source}: {count}")
    else:
        print("No matches found")

def load_proteomexchange_data(xml_file: str = "proteomexchange_all.xml"):
    """Load and clean ProteomeXchange data."""
    xml_path = os.path.join(folder, xml_file)
    
    if not os.path.exists(xml_path):
        print(f"ProteomeXchange XML not found at {xml_path}, downloading...")
        if not download_all_proteomexchange_metadata(xml_path):
            print("Failed to download ProteomeXchange data")
            return None
    else:
        print(f"Using existing ProteomeXchange XML: {xml_path}")
    
    try:
        return create_proteomexchange_dataframe(xml_path)
    except Exception as e:
        print(f"Error processing ProteomeXchange data: {e}")
        return None

def map_pmids(pubmed_path, pride_path, output_path):
    """Legacy function for single mapping - loads PRIDE data each time"""
    pride_df_cleaned = load_and_clean_pride_data(pride_path)
    map_pmids_with_cleaned_pride(pubmed_path, pride_df_cleaned, output_path)

def merge_csv_folder(folder_arg: str, output_name: str = "02_OUTPUT.csv", file_pattern: str = "*pmid_matches.csv") -> Path:
    folder_path = Path(folder_arg)
    if not folder_path.is_dir():
        raise ValueError(f"Not a directory: {folder_path}")
    csv_files = sorted(folder_path.glob(file_pattern))
    if not csv_files:
        raise ValueError(f"No CSV files found in {folder_path}")

    dfs = []
    for p in csv_files:
        try:
            df = pd.read_csv(p, dtype=str, low_memory=False, on_bad_lines="skip")
            df.columns = [c.strip() if isinstance(c, str) else c for c in df.columns]
            dfs.append(df)
        except Exception as e:
            print(f"[WARN] Failed to read {p.name}: {e}", file=sys.stderr)

    if not dfs:
        raise RuntimeError("No readable CSVs found")

    merged = pd.concat(dfs, ignore_index=True, sort=False)

    before = len(merged)
    merged = merged.drop_duplicates(subset=['PMID'], keep="first")
    after = len(merged)
    print(f"Dropped {before - after} duplicate rows by column 'PMID'")

    out_path = Path(output_folder) / output_name
    merged.to_csv(out_path, index=False)
    print(f"Saved merged table to: {out_path}")
    return out_path

if __name__ == "__main__":
    Download()
    Extract()
    
    proteomexchange_df = load_proteomexchange_data("proteomexchange_all.xml")
    
    if proteomexchange_df is None or proteomexchange_df.empty:
        print("⚠️ ProteomeXchange data not available, continuing with PRIDE only")
    else:
        print(f"✓ Loaded {len(proteomexchange_df)} ProteomeXchange records with PMIDs")
    
    if list_of_inputs:
        pubmed_files = sorted(Path(PM_output_folder).glob("*_PMoutput.csv"))
        
        if not pubmed_files:
            print(f"⚠️ No PubMed files (*_PMoutput.csv) found in {PM_output_folder}")
            print("Please run the PubMed download script first (e.g., 01_Literature_table.py)")
            sys.exit(1)
        
        pride_df_cleaned = load_and_clean_pride_data(output_csv_file)
        
        if pride_df_cleaned.empty:
            print("⚠️ No PRIDE data available. Cannot proceed with mapping.")
            sys.exit(1)
        
        print(f"Found {len(pubmed_files)} PubMed files to process.")
        
        for file in pubmed_files:
            print(f"\nProcessing file: {file}")
            map_pmids_output = f"{str(file.name).replace('_PMoutput.csv', '')}_pmid_matches.csv"
            map_pmids_output_path = os.path.join(PM_output_folder, map_pmids_output)
            
            try:
                pubmed_df = pd.read_csv(file)
                if pubmed_df.empty:
                    print(f"No data found in {file}")
                    continue
                
                if proteomexchange_df is not None and not proteomexchange_df.empty:
                    map_pmids_with_proteomexchange_fallback(
                        file, pride_df_cleaned, proteomexchange_df, map_pmids_output_path
                    )
                else:
                    print("ProteomeXchange data not available, using PRIDE only...")
                    map_pmids_with_cleaned_pride(file, pride_df_cleaned, map_pmids_output_path)
                    
            except pd.errors.EmptyDataError:
                print(f"Empty data error for {file}")
            except FileNotFoundError:
                print(f"File not found: {file}")
                continue
                
        merge_csv_folder(folder_arg=output_folder, output_name=main_output)
        
    else:
        main_input_path = os.path.join(output_folder, output_file)
        
        if not os.path.exists(main_input_path):
            print(f"⚠️ Input file not found: {main_input_path}")
            print("Please run the PubMed download script first to generate this file.")
            sys.exit(1)
        
        main_output_path = os.path.join(output_folder, main_output)
        
        try:
            pride_df_cleaned = load_and_clean_pride_data(output_csv_file)
            
            if pride_df_cleaned.empty:
                print("⚠️ No PRIDE data available. Cannot proceed with mapping.")
                sys.exit(1)
            
            if proteomexchange_df is not None and not proteomexchange_df.empty:
                map_pmids_with_proteomexchange_fallback(
                    main_input_path, pride_df_cleaned, proteomexchange_df, main_output_path
                )
            else:
                print("ProteomeXchange data not available, using PRIDE only...")
                map_pmids(main_input_path, output_csv_file, main_output_path)
                
        except pd.errors.EmptyDataError:
            print(f"Empty data error for {main_input_path}")
        except FileNotFoundError:
            print(f"File not found: {main_input_path}")