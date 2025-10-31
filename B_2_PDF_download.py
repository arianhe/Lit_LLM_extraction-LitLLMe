import requests
import pandas as pd
import os
import configparser
from pathlib import Path

def read_config():
    config = configparser.ConfigParser()
    config.read('config.ini')
    return config

if __name__ == "__main__":
    config = read_config()
    output_folder = Path(config.get('Pubmed', 'output_folder'))
    PMID_input = Path(config.get('Input_paths', 'PMID_input'))
    PMID_input = os.path.join(config.get('Pubmed', 'output_folder'), PMID_input)
    PDF_folder = config.get('Input_paths', 'pdf_folder')
    os.makedirs(output_folder / PDF_folder, exist_ok=True)

def get_pdf_url_from_pmid(pmid):
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}%20AND%20SRC:MED&resultType=core&format=json"
    #print(url)
    r = requests.get(url)
    data = r.json()
    try:
        full_texts = data['resultList']['result'][0]['fullTextUrlList']['fullTextUrl']
        for entry in full_texts:
            avail = entry.get('availability', '').lower()
            avail_code = entry.get('availabilityCode', '').lower()
            doc_style = entry.get('documentStyle', '').lower()
            if (
                (avail in ['free', 'open access'] or avail_code == 'oa')
                and doc_style == 'pdf'
            ):
                return entry['url']
    except Exception:
        return None
    return None

def download_pdf(pdf_url, filename):
    try:
        r = requests.get(pdf_url)
        if r.status_code == 200:
            with open(filename, 'wb') as f:
                f.write(r.content)
            return True
        else:
            print(f"HTTP error {r.status_code} for URL: {pdf_url}")
            return False
    except Exception as e:
        print(f"Exception while downloading {pdf_url}: {e}")
        return False

PMID = pd.read_csv(PMID_input)
PMID = PMID.drop_duplicates(subset='PMID', keep='first')
PMID = PMID[PMID['PMID'].notna()]
PMID_list = PMID['PMID'].astype(str).tolist()


i = 0
j = 0
#PMID_list = PMID_list[:2]  # Limit to first 100 PMIDs for testing
for pmid in PMID_list:
    pdf_url = get_pdf_url_from_pmid(pmid)
    if pdf_url:
        save_path = os.path.join(PDF_folder, f"{pmid}.pdf")
        save_path = os.path.join(output_folder, PDF_folder, f"{pmid}.pdf")
        success = download_pdf(pdf_url, save_path)
        if success:
            #print(f"Downloaded {pmid} from {pdf_url}")
            i += 1
        else:
            print(f"Failed to download PDF for PMID {pmid} from {pdf_url}")
            j += 1
    else:
        print(f"No open-access PDF found for PMID {pmid}")
        j += 1
print(f"Total PDFs downloaded: {i}")
print(f"Total PMIDs without open-access PDFs: {j}")

pdf_files = sorted([f.stem for f in Path(output_folder / PDF_folder).glob("*.pdf")])
pmid_in_csv = sorted(PMID_list)
missing_pmids = sorted(set(pmid_in_csv) - set(pdf_files))
pd.Series(missing_pmids).to_csv(f"{output_folder / PDF_folder}/missing_pmids.txt", index=False, header=False)