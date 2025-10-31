#!/usr/bin/env python
# coding: utf-8


import configparser  
import json 
import os
import requests
import pandas as pd
from time import sleep
import unicodedata
import re
from collections import Counter


from pathlib import Path
import Supp_Literatur_table_OpenTargetInt
email, api_key, retmax, output_file, term_list, output_folder, PM_output_folder  = Supp_Literatur_table_OpenTargetInt.read_config()


# Disease search functionality for OpenTargets API
def clean_text(text):
    """Clean text to handle Unicode issues and encoding problems"""
    if text is None:
        return ''

    # Convert to string if not already
    text = str(text)

    # Remove or replace problematic Unicode characters
    # Replace surrogates and other problematic characters
    text = text.encode('utf-8', errors='ignore').decode('utf-8', errors='ignore')

    # Normalize Unicode characters
    try:
        text = unicodedata.normalize('NFKD', text)
    except:
        pass

    # Remove any remaining non-printable characters except common whitespace
    text = ''.join(char for char in text if char.isprintable() or char in [' ', '\t', '\n'])

    # Clean up extra whitespace
    text = ' '.join(text.split())

    return text

# Define GraphQL query to search for diseases
disease_search_query = """
query searchDisease($queryString: String!) {
  search(queryString: $queryString, entityNames: "disease") {
    hits {
      id
      name
      object {
        ... on Disease {
          id
          name
        }
      }
    }
  }
}
"""

# Define GraphQL query to get known drugs for a specific disease
known_drugs_query = """
query KnownDrugsQuery($efoId: String!, $size: Int = 100000) {
  disease(efoId: $efoId) {
    knownDrugs(size: $size) {
      count
      rows {
        phase
        status
        mechanismOfAction
        drug {
          id
          name
          drugType
          maximumClinicalTrialPhase
        }
        target {
          id
          approvedName
          approvedSymbol
        }
      }
    }
  }
}
"""

# Set base URL of GraphQL API endpoint
base_url = "https://api.platform.opentargets.org/api/v4/graphql"

# Search for lung cancer drugs
disease_keywords = term_list
all_results = []
all_drug_names = []  # List to collect all drug names

for disease in disease_keywords:
    print(f"\nSearching for disease: {disease}")

    # Track drug names for this iteration
    iteration_drug_names = []

    try:
        # Step 1: Search for the disease to get its EFO ID
        variables = {"queryString": disease}
        payload = {"query": disease_search_query, "variables": variables}

        r = requests.post(base_url, json=payload, timeout=30)

        if r.status_code != 200:
            print(f"Disease search failed for {disease}")
            print(f"Response: {r.text[:200]}...")
            continue

        search_response = json.loads(r.text)

        # Check for GraphQL errors
        if 'errors' in search_response:
            print(f"GraphQL errors for {disease}: {search_response['errors']}")
            continue

        # Check if we found any disease hits
        hits = search_response.get('data', {}).get('search', {}).get('hits', [])
        if not hits:
            print(f"No disease found for: {disease}")
            continue

        # Get the first (top) hit - this should be the disease EFO ID
        top_hit = hits[0]
        disease_efo_id = top_hit['id']
        disease_name = clean_text(top_hit['name'])  # Clean the disease name

        print(f"Found disease: {disease_name}")

        # Step 2: Get known drugs for this disease
        variables = {"efoId": disease_efo_id, "size": 10000}
        payload = {"query": known_drugs_query, "variables": variables}

        r = requests.post(base_url, json=payload, timeout=30)
        if r.status_code != 200:
            print(f"Known drugs query failed for {disease}")
            print(f"Response: {r.text[:200]}...")
            continue

        drugs_response = json.loads(r.text)

        if 'errors' in drugs_response:
            print(f"GraphQL errors in known drugs query: {drugs_response['errors']}")
            continue

        # Extract known drugs data
        disease_data = drugs_response.get('data', {}).get('disease', {})
        if not disease_data:
            print(f"No disease data found for {disease}")
            continue

        known_drugs = disease_data.get('knownDrugs', {})
        known_drugs_rows = known_drugs.get('rows', [])

        if not known_drugs_rows:
            print(f"No drugs found for {disease}")
            continue

        print(f"Found {len(known_drugs_rows)} drug-target associations for {disease}")
        print(f"Total count: {known_drugs.get('count', 0)}")

        # Process each drug-target association
        for i, drug_row in enumerate(known_drugs_rows):
            drug_data = drug_row.get('drug', {})
            target_data = drug_row.get('target', {})

            if not drug_data:
                print(f"Row {i+1}: Missing drug data")
                continue

            if not target_data:
                print(f"Row {i+1}: Missing target data for drug {drug_data.get('name', 'Unknown')}")
                continue

            # Clean all text fields to handle encoding issues
            drug_name = clean_text(drug_data.get('name', ''))
            target_symbol = clean_text(target_data.get('approvedSymbol', ''))
            target_name = clean_text(target_data.get('approvedName', ''))
            mechanism_of_action = clean_text(drug_row.get('mechanismOfAction', ''))
            status = clean_text(drug_row.get('status', ''))

            # Collect drug name for this iteration
            if drug_name:
                iteration_drug_names.append(drug_name)
                all_drug_names.append(drug_name)  # Add to overall list

            result = {
                'disease_name': disease_name,
                'drug_id': clean_text(drug_data.get('id', '')),
                'drug_name': drug_name,
                'max_clinical_phase': clean_text(str(drug_data.get('maximumClinicalTrialPhase', ''))),
                'status': status,
                'target_id': clean_text(target_data.get('id', '')),
                'target_symbol': target_symbol,
                'target_name': target_name,
                'mechanism_of_action': mechanism_of_action,
            }
            all_results.append(result)

        # Show unique drugs found for this disease
        unique_iteration_drugs = list(set(iteration_drug_names))
        print(f"Successfully processed {len(known_drugs_rows)} drug-target pairs for {disease}")
        print(f"Found {len(unique_iteration_drugs)} unique drugs for {disease}")

    except requests.exceptions.RequestException as e:
        print(f"Network error for {disease}: {e}")
        continue
    except json.JSONDecodeError as e:
        print(f"JSON decode error for {disease}: {e}")
        continue
    except Exception as e:
        print(f"Unexpected error for {disease}: {e}")
        continue

    # Add small delay to be respectful to the API
    sleep(1)

# Count drug name frequencies and sort by count
drug_counter = Counter([clean_text(name) for name in all_drug_names if clean_text(name)])
unique_drug_names = list(drug_counter.keys())

print(f"\nProcessing complete!")
print(f"Total drug-target associations collected: {len(all_results)}")
print(f"Total drug names collected (with duplicates): {len(all_drug_names)}")
print(f"Unique drug names (deduplicated): {len(unique_drug_names)}")

# Convert results to DataFrame and display results (OUTSIDE the for loop)
if len(all_results) > 0:
    df_results = pd.DataFrame(all_results)

    # Save all results to CSV with proper encoding
    output_file = "01_all_drugs_opentargets.csv"
    try:
        df_results.to_csv(Path(PM_output_folder) / output_file, index=False, encoding='utf-8')
        print(f"\nAll results saved to: {PM_output_folder}{output_file}")
    except UnicodeEncodeError as e:
        print(f"Encoding error when saving CSV, trying with error handling: {e}")
        df_results.to_csv(Path(PM_output_folder) / output_file, index=False, encoding='utf-8', errors='ignore')
        print(f"All results saved to: {PM_output_folder}{output_file} (with encoding error handling)")

    # Show kinase-related drugs if any
    kinase_related = df_results[df_results['target_name'].str.contains('kinase', case=False, na=False)]

    if len(kinase_related) > 0:
        print(f"\n🧬 Found {len(kinase_related)} kinase-related entries:")
        print(f"🧬 Unique kinase-related drugs: {len(kinase_related['drug_id'].unique())}")
        print(f"🧬 Unique kinase targets: {kinase_related['target_symbol'].nunique()}")

        # Save kinase-related results separately
        try:
            kinase_related.to_csv(Path(PM_output_folder) / "01_lung_cancer_kinase_related_drugs_opentargets.csv", index=False, encoding='utf-8')
            print(f"Kinase-related results saved to: 01_lung_cancer_kinase_related_drugs_opentargets.csv")
        except UnicodeEncodeError:
            kinase_related.to_csv(Path(PM_output_folder) / "01_lung_cancer_kinase_related_drugs_opentargets.csv", index=False, encoding='utf-8', errors='ignore')
            print(f"Kinase-related results saved to: 01_lung_cancer_kinase_related_drugs_opentargets.csv (with encoding error handling)")

        kinase_list = kinase_related['drug_name'].unique().tolist()
        # Save to parquet for LLM
        kinase_related.to_parquet(Path(PM_output_folder) / "01_lung_cancer_kinase_related_drugs_opentargets.parquet", index=False)

        # Count kinase drug frequencies specifically
        kinase_drug_counter = Counter(kinase_related['drug_name'].tolist())

        # Save kinase drugs with counts
        kinase_drug_df = pd.DataFrame([
            {'drug_name': drug, 'count': count} 
            for drug, count in kinase_drug_counter.most_common()
        ])
        try:
            kinase_drug_df.to_csv("kinase_drugs_with_counts.csv", index=False, encoding='utf-8')
            print(f"Kinase drugs with counts saved to: kinase_drugs_with_counts.csv")
        except UnicodeEncodeError:
            kinase_drug_df.to_csv("kinase_drugs_with_counts.csv", index=False, encoding='utf-8', errors='ignore')
            print(f"Kinase drugs with counts saved to: kinase_drugs_with_counts.csv (with encoding error handling)")

        # Display kinase drugs by frequency
        print(f"\nKinase-related drugs sorted by frequency:")
        for i, (drug, count) in enumerate(kinase_drug_counter.most_common(), 1):
            print(f"  {i:2d}. {drug} (appears {count} times)")

    else:
        print("\nNo kinase-related drugs found")
        kinase_list = []  # Define empty list if no kinase drugs found

    # Save the unique drug names list with counts (all drugs)
    unique_drug_df = pd.DataFrame([
        {'drug_name': drug, 'count': count} 
        for drug, count in drug_counter.most_common()
    ])
    
else:
    print("No results found for any diseases")
    kinase_list = []  # Define empty list if no results found



email, api_key, retmax, output_file, term_list, output_folder, PM_output_folder  = Supp_Literatur_table_OpenTargetInt.read_config()
Supp_Literatur_table_OpenTargetInt.Literature_search(kinase_list, email, api_key, output_folder, output_file, PM_output_folder, retmax)
Supp_Literatur_table_OpenTargetInt.merge_csv_folder(PM_output_folder=PM_output_folder, output_name=output_file, output_folder=output_folder)
