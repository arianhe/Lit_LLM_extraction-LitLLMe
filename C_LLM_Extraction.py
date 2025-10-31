import os
import re
import json
import time
import difflib
import traceback
import pandas as pd
from tqdm import tqdm #Progressbar
from llama_cpp import Llama
import ast #Postprocessing
import configparser
from pathlib import Path
import fitz  # PyMuPDF for PDF processing
from B_1_ProtXchange import merge_csv_folder 
print("Packages loaded successfully.")



def read_config():
    config = configparser.ConfigParser()
    config.read('config.ini')
    return config

if __name__ == "__main__":
    config = read_config()
    # Access variables directly from any section
    tissue_list = json.loads(config.get('Pubmed', 'term_list'))
    list_output_folder = config.get('Pubmed', 'list_output_folder')
    output_folder = config.get('Pubmed', 'output_folder')
    input_folder = os.path.join(output_folder, list_output_folder)
    os.makedirs(input_folder, exist_ok=True)
    PM_output_folder = input_folder
    list_of_inputs = config.getboolean('Input_paths', 'List_of_Inputs')
    file_ending = config.get('Input_paths', 'File_ending')
    main_input_LLM = config.get('Input_paths', 'Main_input_LLM')
    main_input_LLM = os.path.join(output_folder, main_input_LLM)
    kinase_input = config.get('Input_paths', 'Kinases_input')
    cell_input = config.get('Input_paths', 'Cell_input')
    model_path = config.get('LLM', 'model_path')
    device_id = int(config.get('LLM', 'device_id'))
    n_ctx = int(config.get('LLM', 'n_ctx'))
    n_threads = int(config.get('LLM', 'n_threads'))
    batch_size = int(config.get('LLM', 'batch_size'))
    output_file = config.get('Output', 'output_file')
    PDF_extraction = config.getboolean('Input_paths', 'PDF_extraction')
    PDF_output_file = config.get('Input_paths', 'PDF_output_file')
    paperfolder = Path(os.path.join(output_folder, 'Papers'))

    


print("Configuration loaded successfully.")

if device_id not in [-1, 0]: raise ValueError("Only device_ids 0 and -1 are allowed")

# === Model loading ===
def load_llm(device_id: int, model_path: str):
    if device_id != 0:
        print(f"Using device 0 (TITAN Xp) is recommended. If device {device_id} is not a GPU, please reconsider....")
        time.sleep(5)

    #os.environ["CUDA_VISIBLE_DEVICES"] = str(device_id)
    n_gpu_layers = 100  # TITAN Xp full offload
    if device_id == 0:
        try:
            llm_model = Llama(
                model_path=model_path,
                device=device_id,
                n_ctx=n_ctx,
                n_threads=n_threads,
                n_gpu_layers=n_gpu_layers,
                verbose=False
            )
            print(f"Model successfully loaded on GPU {device_id} (TITAN Xp)")
            return llm_model

        except Exception as e:
            print(f"Model loading on Device {device_id} failed. Try devide_id = -1")
            traceback.print_exc()
    if device_id == -1:
        try:
            llm_model = Llama(
                model_path=model_path,
                device=device_id,  # -1 for CPU
                n_ctx=4096,
                n_threads=n_threads,
                verbose=False
            )
            print("Model successfully loaded on CPU.")
            return llm_model
        except Exception as cpu_e:
            print("Model loading on CPU also failed.")
            traceback.print_exc()
            raise cpu_e

def fix_inner_quotes_in_lists(text):
    def repl(match):
        items = match.group(1).split(',')
        cleaned_items = []
        for item in items:
            item = item.strip()
            # Only escape double quotes inside the item, not at the boundaries
            if item.startswith('"') and item.endswith('"'):
                inner = item[1:-1].replace('"', '\\"')
                cleaned_items.append(f'"{inner}"')
            else:
                cleaned_items.append(item.replace('"', '\\"'))
        return '[' + ', '.join(cleaned_items) + ']'
    return re.sub(r'\[([^\[\]]+?)\]', repl, text)


# === JSON Cleaning Utility ===
def clean_json_text(text: str) -> str:
    text = text.encode("utf-8", errors="ignore").decode("utf-8", errors="ignore")
    text = re.sub(r'[\x00-\x1f]+', '', text)
    text = re.sub(r',\s*([\]}])', r'\1', text)
    # Replace double quotes inside lists with single quotes
    text = fix_inner_quotes_in_lists(text)
    # Fix quoted lists like "[a, b, c]" -> ["a", "b", "c"]
    def fix_stringified_list(match):
        items = [item.strip().strip('"') for item in match.group(1).split(',')]
        quoted_items = [json.dumps(i) for i in items if i]
        return "[" + ", ".join(quoted_items) + "]"
    text = re.sub(r'"\[([^\[\]]+?)\]"', fix_stringified_list, text)

    # Fix dict-like sets {a, b, c} → ["a", "b", "c"]
    def fix_set_like(match):
        items = [i.strip().strip('"') for i in match.group(1).split(',')]
        quoted_items = [json.dumps(i) for i in items if i]
        return "[" + ", ".join(quoted_items) + "]"
    text = re.sub(r'\{([^{}:]+?)\}', fix_set_like, text)

    # Quote keys like Drug: → "Drug":
    def quote_keys(m):
        key = m.group(1)
        if key.startswith('"') and key.endswith('"'):
            return m.group(0)
        return f'"{key}":'
    text = re.sub(r'(\w+)\s*:', quote_keys, text)

    # Fix stray backslashes
    text = re.sub(r'\\(?!["\\/bfnrtu])', r'\\\\', text)

    # Final cleanup and ensure it ends with a full brace
    text = re.sub(r',\s*([\]}])', r'\1', text)
    if text.count('[') > text.count(']'):
        text += ']'
    if text.count('{') > text.count('}'):
        text += '}'

    return text.strip()

def try_partial_json_parse(text: str) -> dict:
    start = text.find('{')
    end = text.rfind('}')
    if start != -1 and end != -1 and end > start:
        candidate = text[start:end+1]
        try:
            parsed = json.loads(candidate)
        except Exception:
            cleaned = clean_json_text(candidate)
            try:
                parsed = json.loads(cleaned)
            except Exception:
                return {}

        # Postprocess: dedupe + limit list length
        for field in ["Drug", "Cell Type", "Kinase"]:
            if field in parsed and isinstance(parsed[field], list):
                unique_items = list(set(parsed[field]))
                parsed[field] = unique_items[:20]  # cap to first 20
        return parsed
    return {}


def extract_field(text, field_name):
    pattern = rf'"{re.escape(field_name)}"\s*:\s*("(.*?)"|(\[.*?\])|(\w+))'
    match = re.search(pattern, text, re.DOTALL)
    if not match:
        return ""
    val = match.group(2) or match.group(3) or match.group(4) or ""
    return val.strip('"').strip() if isinstance(val, str) else val

def detect_cell_lines(abstract, known_lines):
    words = abstract.split()
    found = []
    for cl in known_lines:
        if cl in abstract or difflib.get_close_matches(cl, words, n=1, cutoff=0.9):
            found.append(cl)
    return found

def extract_metadata(abstract: str, kinase_inhibitors: list, cell_lines: list) -> dict:
    global llm

    detected_cell_lines = detect_cell_lines(abstract, cell_lines)
    drugs_list = ', '.join(sorted(kinase_inhibitors))
    kinase_list = ','.join(sorted(kinases))
    

    with open("Prompt.json", "r", encoding="utf-8") as f:
        loaded_messages = json.load(f)


    system_content = loaded_messages[0]["content"]
    user_content = loaded_messages[1]["content"].format(drugs_list, kinase_list, detected_cell_lines, abstract)

    messages = [
        {"role": "system", "content": system_content},
        {"role": "user", "content": user_content}
    ]

    try:
        output = llm.create_chat_completion(messages=messages, max_tokens=1024)
    except Exception as e:
        print(f"[LLM ERROR] Failed to call model: {e}")
        return {}

    if not output or "choices" not in output:
        print(f"[LLM ERROR] Output missing 'choices': {repr(output)}")
        return {}

    try:
        raw_text = output["choices"][0]["message"]["content"]
    except Exception as e:
        print(f"[LLM ERROR] Unexpected format in response: {e}\n→ Output: {repr(output)}")
        return {}

    for attempt in range(3):
        try:
            cleaned = clean_json_text(raw_text)
            metadata = json.loads(cleaned)

            if not isinstance(metadata, dict):
                raise ValueError("Parsed JSON is not a dictionary")

            # Normalize fields
            metadata["Drug"] = str(metadata.get("Drug", "")).strip()

            val = metadata.get("Cell Type", "")
            if val and not isinstance(val, list):
                metadata["Cell Type"] = [val]
            elif not val:
                metadata["Cell Type"] = []

            return metadata

        except Exception as e:
            print(f"[Parse Attempt Failed #{attempt + 1}] Error: {e}")
            print(f"[Cleaned JSON String]: {cleaned}")
            continue

    # Fallback: basic regex-based extraction
    print(f"[Parse Warning] Could not fully decode JSON, falling back.")

    fallback_metadata = {
    "Drug": extract_field(raw_text, "Drug") or "",
    "Cell Type": parse_list_col_safe(extract_field(raw_text, "Cell Type")) or [],
    "Kinase": parse_list_col_safe(extract_field(raw_text, "Kinase")) or []
}

    print(f"[Fallback Extracted Metadata]: {fallback_metadata}")
    return fallback_metadata

def process_batch(args):
    batch_df, kinase_inhibitors, cell_lines = args
    results = []

    for idx, row in batch_df.iterrows():
        if PDF_extraction == True:
            abstract = row["Methods_Text"]
        else:  
            abstract = row["ExtractionColumn"]
        metadata = {}

        for attempt in range(3):
            try:
                metadata = extract_metadata(abstract, kinase_inhibitors, cell_lines)
                break
            except RuntimeError as e:
                if "cudaMalloc" in str(e) or "out of memory" in str(e).lower():
                    print(f"[OOM] Skipping abstract at index {idx}.")
                    break
                time.sleep(1 + attempt)
            except Exception as e:
                if attempt == 2:
                    print(f"[ERROR] Unexpected failure at index {idx}: {e}")
                    print(f"[DEBUG] Abstract snippet: {abstract[:200].replace('\n', ' ')}...")
                    break

        results.append({
            'index': idx,
            'Drug': metadata.get("Drug", ""),
            'Cell Type': metadata.get("Cell Type", []),
            'Kinase': metadata.get("Kinase", [])
        })

    return results

def extract_metadata_serial(df, kinase_inhibitors, cell_lines, batch_size=10):
    all_batches = [df.iloc[i:i+batch_size] for i in range(0, len(df), batch_size)]
    results = []

    print(f"Processing {len(all_batches)} batches serially on device {device_id}")

    for batch in tqdm(all_batches, desc="Extracting metadata"):
        batch_result = process_batch((batch, kinase_inhibitors, cell_lines))
        results.append(batch_result)

    flat_records = [record for batch in results for record in batch]
    return pd.DataFrame(flat_records).set_index("index")

#Post process cleaning
def parse_list_col_safe(x):
    # If it’s already a list, return as is
    if isinstance(x, list):
        return x
    # If it’s a string looking like a list, parse it safely
    if isinstance(x, str):
        x = x.strip()
        if x.startswith('[') and x.endswith(']'):
            try:
                return ast.literal_eval(x)
            except Exception:
                # fallback if parsing fails
                return []
        # Handle empty strings or commas or weird formats
        if x in ['', ',']:
            return []
        # If it's a single item string, put it in a list
        return [x]
    # For NaNs or anything else, return empty list to avoid crashing
    return []

def extract_methods_section(pdf_path):
    doc = fitz.open(pdf_path)
    methods_text = ""
    in_methods = False
    methods_found = False

    heading_pattern = re.compile(
    r"^[^\w]*"  # Any non-word characters (including ■, spaces, etc.)
    r"(?:\d+(\.\d+)*\s*[\.\)]?\s*)?"  # Optional numbering
    r"(materials?\s*(and|&)?\s*)?methods\b"
    r"|materials\b"  # <-- Add | to make it an alternative
    r"|star\s*[^\w]*methods\b"
    r"|experimental\s*section\b"
    r"|online\s*methods\b",
    re.IGNORECASE
)
    
    stop_headings = [
        "results", "discussion", "conclusion", "references", "acknowledgements", "introduction", "SUPPLEMENTARY"
    ]
    
    stop_pattern = re.compile(
        r"^(?:\d+(\.\d+)*\s*[\.\)]?\s*)?(" + "|".join(stop_headings) + r")(\b|:| and|/|s)?", re.IGNORECASE
    )

    for page in doc:
        blocks = page.get_text("blocks")
        for b in blocks:
            text = b[4].strip()
            # Find methods section
            if not methods_found and heading_pattern.search(text):
                in_methods = True
                methods_found = True
                match = heading_pattern.search(text)
                start_pos = match.start()
                methods_text += "\n" + text[start_pos:]
                continue
                
            # Check for end of methods section
            if in_methods and stop_pattern.match(text):
                in_methods = False
                
            # Add text if we're in the methods section
            if in_methods:
                methods_text += "\n" + text

    # Print status ONCE after processing all pages
    if methods_found:
        print(f"✅ Methods section found in {pdf_path}")

    return methods_text.strip()

print("Methods loaded successfully")

llm = load_llm(device_id=device_id, model_path=model_path) 


#Helpful ressources
try:
    kinase_table = pd.read_parquet(Path(PM_output_folder) / "01_lung_cancer_kinase_related_drugs_opentargets.parquet")
    kinase_inhibitors = kinase_table['drug_name'].dropna().unique().tolist()
    kinases = kinase_table['target_name'].dropna().unique().tolist()
    print(f"Kinase table loaded from parquet {PM_output_folder}01_lung_cancer_kinase_related_drugs_opentargets.parquet.")
except Exception as e:
    print(f"Error reading kinase table from parquet {PM_output_folder}01_lung_cancer_kinase_related_drugs_opentargets.parquet: {e}")
    kinase_table = pd.read_csv(kinase_input, sep=';')  # Fallback to CSV
    kinase_inhibitors = kinase_table['Inhibitor'][:].to_list()
    kinases = kinase_table['Target'][:].to_list()
    print(f"Kinase table loaded from CSV {kinase_input}.")
cell_lines_base = pd.read_csv(cell_input)
cell_lines = cell_lines_base['name'][:].to_list()

#Structure of future dataframe

start = time.perf_counter()
if list_of_inputs == False:
    df = pd.read_csv(main_input_LLM)
    print(f"Input file {main_input_LLM} read successfully with {len(df)} entries.")
    #Delete duplicates decided by pmid
    df = df.drop_duplicates(subset='PMID', keep='first')
    
    # Remove rows with NaN or empty/whitespace-only abstracts
    df = df[df['Abstract'].fillna('').str.strip() != '']
    Test = df.copy()
    Test['Drug'] = ''
    Test['Cell Type'] = ''
    Test['Kinase'] = ''
    Test['ExtractionColumn'] = Test['Abstract'].fillna('') + ' ' + Test['sampleProcessingProtocol'].fillna('')


    print(f"Processing {len(Test)} abstracts with {len(kinase_inhibitors)} kinase inhibitors and {len(cell_lines)} cell lines.")
    df_results = extract_metadata_serial(
        df=Test,
        kinase_inhibitors=kinase_inhibitors,
        cell_lines=cell_lines,
        batch_size=batch_size  
)
    for col in ['Drug','Cell Type']:
        df_results[col] = df_results[col].apply(parse_list_col_safe)

    Test.update(df_results)
    Test.drop(columns=['ExtractionColumn'], inplace=True)
    Test.to_csv(f"{output_folder}{output_file}", index=False)


elif PDF_extraction == False and list_of_inputs == True:
    for tissue in tissue_list:
        tissue = f"{tissue.replace(' ', '_')}"
        try:
            df = pd.read_csv(f"{input_folder}{tissue}{file_ending}")
            if df.empty == 1:
                print(f"⚠️ ValueError for {tissue}: DataFrame is empty")
                continue
            #Delete duplicates decided by pmid
            df = df[df['PMID_clean'].notna()]
            df = df.drop_duplicates(subset=['PMID_clean'])
            # Remove rows with NaN or empty/whitespace-only abstracts
            df = df[df['Abstract'].fillna('').str.strip() != '']
            Test = df.copy()
            Test['Drug'] = ''
            Test['Cell Type'] = ''
            Test['Kinase'] = ''
            Test['ExtractionColumn'] = Test['Abstract'].fillna('') + ' ' + Test['sampleProcessingProtocol'].fillna('')
            print(f"Processing file: {tissue}")
            print(f"Processing {len(Test)} abstracts with {len(kinase_inhibitors)} kinase inhibitors and {len(cell_lines)} cell lines.")
            df_results = extract_metadata_serial(
                df=Test,
                kinase_inhibitors=kinase_inhibitors,
                cell_lines=cell_lines,
                batch_size=batch_size  
            )
            for col in ['Drug','Cell Type']:
                df_results[col] = df_results[col].apply(parse_list_col_safe)
            output_file = f"{input_folder}{tissue}_extracted.csv"
            Test.update(df_results)
            Test.drop(columns=['ExtractionColumn'], inplace=True)
            Test.to_csv(output_file, index=False)
        except pd.errors.EmptyDataError:
            print(f"⚠️ No data found for {tissue}. Skipping.")
            continue
        except FileNotFoundError:
            print(f"⚠️ File not found for {tissue}. Skipping.")
            continue
elif PDF_extraction == True:
    data = []
    for paper in paperfolder.glob("*.pdf"):
    #paperlist = ['Outputs/Papers/20190765.pdf', 'Outputs/Papers/21368869.pdf']  # Example list of PDF files
    #for paper in paperlist:
        file_base = os.path.splitext(os.path.basename(paper))[0]
        try:
            methods_text = extract_methods_section(paper)
            if not methods_text.strip():
                print(f"⚠️ No methods section found in {paper}. Skipping.")
        except Exception as e:
            print(f"⚠️ Error processing {paper}: {e}")
            methods_text = ""
        data.append({'Filename': file_base, 'Methods_Text': methods_text})
    methods_df = pd.DataFrame(data)
    df = methods_df    
    if df.empty == 1:
        print(f"⚠️ ValueError: DataFrame is empty")
    else:
        #Delete duplicates decided by filename
        df = df.drop_duplicates(subset='Filename', keep='first')
        # Remove rows with NaN or empty/whitespace-only abstracts
        df = df[df['Methods_Text'].fillna('').str.strip() != '']
        df['Drug'] = ''
        df['Cell Type'] = ''
        df['Kinase'] = ''

        df_results = extract_metadata_serial(
            df=df,
            kinase_inhibitors=kinase_inhibitors,
            cell_lines=cell_lines,
            batch_size=batch_size  
        )
        for col in ['Drug','Cell Type', 'Kinase']:
            df_results[col] = df_results[col].apply(parse_list_col_safe)
        os.makedirs(paperfolder, exist_ok=True)
        PDF_output_path = os.path.join(paperfolder, PDF_output_file)
        df.update(df_results)    
        df.to_csv(PDF_output_path, index=False)
        print(f"✅ PDF extraction completed. Results saved to {PDF_output_path}")

end = time.perf_counter()
print(f"✅ Finished in {(end - start)/60:.2f} minutes")