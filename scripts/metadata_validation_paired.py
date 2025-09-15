########################################## Loading libraries and user arguments ##########################################

import pandas as pd
import os
import sys
import logging
import csv

### TO DO:
# - [OK] Verify the presence of a subtitle row
# - [OK] Check for '#q2:types' in the first row, and handle sample-id column case variations
# - [OK] Detect and validate 'categorical' and 'numeric' subtitle types
# - [OK] Strip whitespace from all cells
# - [OK] Replace '/' characters with '-' in sample names
# - [OK] Ensure numeric columns only contain int or float values
# - [OK] Verify that sequence names match sample-ids in manifest.tsv
# - [OK] Confirm matching sample counts between Manifest and Metadata files
# - [OK] Add support for paired-end data:
#        - If param_paired = True, sample-ids in Metadata should correspond to directories
#        - If param_paired = False, sample-ids should correspond to files
# - [TODO] Improve use of the logging library for better control of messages, warnings, and errors
# - [OK] When paired-end: integrate manifest file (from fastqc_validation_paired.py) with Metadata
#        - Use directory names as sample-ids
#        - Use directory paths as absolute-filepaths
#        - Remove duplicate entries

########################################## Loading libraries and user arguments ##########################################

import pandas as pd
import os
import sys
import logging
import csv
import argparse

### TO DO:
# - [OK] Verify the presence of a subtitle row
# - [OK] Check for '#q2:types' in the first row, and handle sample-id column case variations
# - [OK] Detect and validate 'categorical' and 'numeric' subtitle types
# - [OK] Strip whitespace from all cells
# - [OK] Replace '/' characters with '-' in sample names
# - [OK] Ensure numeric columns only contain int or float values
# - [OK] Verify that sequence names match sample-ids in manifest.tsv
# - [OK] Confirm matching sample counts between Manifest and Metadata files
# - [OK] Add support for paired-end data:
#        - If param_paired = True, sample-ids in Metadata should correspond to directories
#        - If param_paired = False, sample-ids should correspond to files
# - [TODO] Improve use of the logging library for better control of messages, warnings, and errors
# - [OK] When paired-end: integrate manifest file (from fastqc_validation_paired.py) with Metadata
#        - Use directory names as sample-ids
#        - Use directory paths as absolute-filepaths
#        - Remove duplicate entries

# Set up argument parser
parser = argparse.ArgumentParser(description='Process metadata and samples.')
parser.add_argument('--test', action='store_true', help='Run in development mode with hardcoded paths.')

# Define command-line arguments
parser.add_argument('--input_metadata_path', type=str, help='Path to the input metadata CSV file.')
parser.add_argument('--input_samples_directory', type=str, help='Directory containing input samples.')
parser.add_argument('--input_manifest_path', type=str, help='Path to the manifest TSV file.')
parser.add_argument('--param_user_separator', type=str, default='\t', help='User-specified separator for files (default: \\t).')
parser.add_argument('--param_sample_identifier', type=str, default='LABID', help='Column name to use as sample identifier (default: LABID).')
parser.add_argument('--param_fill_nan_values', type=str, default='median', help="Method to fill NaN values: 'median', 'mean', or a float (default: median).")
parser.add_argument('--output_validated_metadata_tsv', type=str, help='Output path for validated metadata TSV.')
parser.add_argument('--output_report_metadata_html', type=str, help='Output path for the report HTML.')
parser.add_argument('--output_logs', type=str, help='Output path for logs.')
parser.add_argument('--param_paired', type=lambda x: x.lower() == 'true', default=True, help='Whether the data is paired-end: true or false (default: true).')

args = parser.parse_args()

test=False
if args.test or test:
    print('Development mode activated!')
    # Hardcoded paths and parameters for development mode
    input_metadata_path = '/home/jpereira/shotgun_snake/Metadata/metadata_metag_2025-4-10.csv'
    input_samples_directory = '/home/pperaza/NOVOGENE/MG/20240711_Respaldo_Pablo/01.RawData/'
    input_manifest_path = '/home/jpereira/shotgun_snake/Results_test/Data/Fastqc/manifest.tsv'
    param_user_separator = '\t'
    param_sample_identifier = 'LABID'
    param_fill_nan_values = 'median'
    param_paired = True
    output_validated_metadata_tsv = '/home/jpereira/shotgun_snake/Results_test/Data/Metadata/internal_metadata.tsv'
    output_report_metadata_html = '/home/jpereira/shotgun_snake/Results_test/Report/Metadata/internal_metadata.html'
    output_logs = '/home/jpereira/shotgun_snake/Results_test/Data/Metadata/logs.txt'
else:
    # Use command-line provided arguments
    input_metadata_path = args.input_metadata_path
    input_samples_directory = args.input_samples_directory
    input_manifest_path = args.input_manifest_path
    param_user_separator = args.param_user_separator
    param_sample_identifier = args.param_sample_identifier
    param_fill_nan_values = args.param_fill_nan_values
    param_paired = args.param_paired
    output_validated_metadata_tsv = args.output_validated_metadata_tsv
    output_report_metadata_html = args.output_report_metadata_html
    output_logs = args.output_logs


# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(output_logs, mode='w'),
        logging.StreamHandler(sys.stdout)
    ]
)

errors = False

###### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(output_logs, mode='w'),
        logging.StreamHandler(sys.stdout)
    ]
)

errors = False

###### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

############################ Defining metadata file separator, loading data, filling NaN values ############################

def detect_file_type(file_path):
    with open(file_path, 'r') as file:
        sample = file.read(4096)  # Read a sample of the file
    sniffer = csv.Sniffer()
    try:
        dialect = sniffer.sniff(sample)
        return dialect.delimiter
    except csv.Error:
        logging.error("Could not detect delimiter automatically. The metadata file might have an unrecognized format.")
        return None

dialect_delimiter = detect_file_type(input_metadata_path)

# Determine final delimiter
if param_user_separator in ['', None]:
    separator = dialect_delimiter
else:
    separator = param_user_separator

# Raise an error if the final delimiter is '' or None
if separator in ['', None]:
    logging.error("Error: Values such as '' or 'None' cannot be used as a delimiter for the Metadata file")
    exit(1)

# Warn if user-provided separator does not match the auto-detected one
if separator != dialect_delimiter:
    logging.warning(
        f"Warning: Column delimiter provided by the user '{separator}' "
        f"does not match the one automatically detected '{dialect_delimiter}'"
    )

# Read the metadata
metadata_df = pd.read_csv(input_metadata_path, sep=separator)

# Strip leading and trailing whitespace from all string cells
metadata_df = metadata_df.map(lambda x: x.strip() if isinstance(x, str) else x)

# Replace '/' characters with '-' characters in all string cells
metadata_df = metadata_df.map(lambda x: x.replace('/', '-') if isinstance(x, str) else x)

# Fill NaN values based on user parameter
if param_fill_nan_values == 'median':
    metadata_df = metadata_df.fillna(metadata_df.median(numeric_only=True))
elif param_fill_nan_values == 'mean':
    metadata_df = metadata_df.fillna(metadata_df.mean(numeric_only=True))
elif isinstance(param_fill_nan_values, (int, float)):
    metadata_df = metadata_df.fillna(param_fill_nan_values)
else:
    logging.error(
        "Error: Unknown value for argument 'fill_nan_values'. "
        "Allowed values are: 'median', 'mean', float, or int."
    )
    exit(1)

###### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

############### Checking the presence of subtitles, making a standard metadata file with QIIME2-like format ################

# Check for comments and QIIME2 subtitles
first_col = metadata_df.iloc[:,0]

# Check for comments
comment_bool_list = first_col.astype(str).str.contains(r'^#')
comment_index_list = comment_bool_list[comment_bool_list].index.to_list()

# Check presence of comment rows
if comment_index_list:
    logging.info(f"Comments detected in rows: {comment_index_list}")
    comments = True
else:
    logging.info("No comment rows detected")
    comments = False

# Check presence of subtitle rows
if 0 in comment_index_list and first_col[0] == '#q2:types':
    logging.info("Subtitle row detected")
    subtitle = True
else:
    logging.info("No subtitle row detected")
    subtitle = False

# Check for 'categorical' and 'numeric' categories in subtitle row
if subtitle:
    subtitle_categories = set(metadata_df.iloc[0, :])
    unrecognized_categories = subtitle_categories - set(['#q2:types', 'numeric', 'categorical'])
    if unrecognized_categories:
        logging.error(f"Error: The following subtitle categories are not recognized: {unrecognized_categories}")
        logging.error("Only '#q2:types', 'numeric', and 'categorical' can be used as categories in the subtitle row")
        errors = True

# Infer data types of the metadata_df:
if not subtitle:
    if param_sample_identifier in ['', None]:
        logging.error("Error: No subtitle row was found and no sample identifier column was provided.")
        logging.error("Please define a subtitle row or provide the column name with the sample files names.")
        errors = True
    elif param_sample_identifier not in metadata_df.columns:
        logging.error("Error: The sample identifier provided does not match any of the column headers in the metadata.")
        logging.error(f"Column headers: {metadata_df.columns}")
        errors = True

# Validation of 'numeric' columns if they are present
if subtitle:
    for col in metadata_df.columns:
        if metadata_df[col][0] == 'numeric':
            true_nums = pd.to_numeric(metadata_df[col][1:], errors='coerce').notna()
            
            # Check if any value in true_nums is False
            if not true_nums.all():
                logging.error(f"Error: Column {col} with subtitle 'numeric' seems to have cells with string values")
                # Identify rows where true_nums is False
                offending_rows = true_nums[~true_nums].index
                logging.error(f"Rows with offending strings: {offending_rows.tolist()}")
                errors = True

if errors:
    logging.error("Errors detected when trying to format the metadata file. Exiting script.")
    exit(1)

# Generalizing metadata file format
if subtitle:
    metadata_df = metadata_df.rename(columns={metadata_df.columns[0]: 'sample-id'})
else:
    # Rename the user-provided column to 'sample-id'
    metadata_df = metadata_df.rename(columns={param_sample_identifier: 'sample-id'})

# Move 'sample-id' column to the first position
cols = ['sample-id'] + [col for col in metadata_df.columns if col != 'sample-id']
metadata_df = metadata_df[cols]

# Adding a row with subtitles if not present
if not subtitle:    
    logging.info("Adding subtitles to the first row:")
    subtitle_row = []
    for col in metadata_df.columns:
        if col == 'sample-id':
            subtitle_row.append('#q2:types')
        else:
            all_nums = pd.to_numeric(metadata_df[col], errors='coerce').notna().all()
            if all_nums:
                subtitle_row.append('numeric')
                metadata_df[col] = pd.to_numeric(metadata_df[col], errors='coerce')
                logging.info(f"Adding 'numeric' subtitle to column: {col}")
            else:
                subtitle_row.append('categorical')
                logging.info(f"Adding 'categorical' subtitle to column: {col}")

subtitle_df = pd.DataFrame([subtitle_row], columns=metadata_df.columns)
metadata_df = pd.concat([subtitle_df, metadata_df], axis=0)

###### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

############################################## Checking sample-id sequences ##############################################

# Check if sample-ids from the metadata are present as directories or files in input_samples_directory
sample_error = False

if param_paired:
    logging.warning('param_paired = True: Checking if each sample-id corresponds to a directory in the input samples directory.')
    logging.warning('Manifest file checking will be skipped for paired-end data.')
else:
    logging.warning('param_paired = False: Checking if each sample-id corresponds to a file in the input samples directory.')

for sample in metadata_df['sample-id'].iloc[1:]:
    sample_path = os.path.join(input_samples_directory, sample)
    
    if param_paired:
        if not os.path.isdir(sample_path):
            logging.error(f"Error: Sample-id '{sample}' does not correspond to a directory.")
            logging.error(f"Expected directory path: '{sample_path}'")
            sample_error = True
    else:
        if not os.path.isfile(sample_path):
            logging.error(f"Error: Sample-id '{sample}' does not correspond to a file.")
            logging.error(f"Expected file path: '{sample_path}'")
            sample_error = True

if sample_error:
    exit(1)

# Read manifest and preprocess for paired-end data
manifest_df = pd.read_csv(input_manifest_path, sep='\t', header=0)

if param_paired:
    manifest_df['absolute-filepath'] = manifest_df['absolute-filepath'].apply(os.path.dirname)
    manifest_df['sample-id'] = manifest_df['absolute-filepath'].apply(os.path.basename)
    manifest_df = manifest_df.drop_duplicates()

# Create subtitle row and prepend it
subtitle_row = pd.DataFrame([['#q2:types'] + ['categorical'] * (manifest_df.shape[1] - 1)],
                            columns=manifest_df.columns)
manifest_df = pd.concat([subtitle_row, manifest_df], ignore_index=True)

# Compare sample presence in metadata and manifest
missing_in_metadata = set(manifest_df['sample-id']) - set(metadata_df['sample-id'])
missing_in_manifest = set(metadata_df['sample-id']) - set(manifest_df['sample-id'])

if missing_in_metadata:
    logging.warning(f"Warning: These samples are present in the manifest but missing from the metadata:\n{missing_in_metadata}")
if missing_in_manifest:
    logging.warning(f"Warning: These samples are present in the metadata but not found in the manifest:\n{missing_in_manifest}")

# Merge to retain only shared samples
metadata_df = pd.merge(manifest_df, metadata_df, on='sample-id', how='inner')

###### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

###################################### Saving validated metadata in tsv and html files  #####################################
logging.info(f'Saving metadata file in: {output_validated_metadata_tsv}')
os.makedirs(os.path.dirname(output_validated_metadata_tsv), exist_ok=True)
metadata_df.to_csv(output_validated_metadata_tsv, sep='\t', index=False)

if not metadata_df.empty:
    html = metadata_df.to_html(escape=False, classes='table table-striped table-hover')
else:
    html = "<p>No data available in the report.</p>"

logging.info(f'Saving metadata file in html format in: {output_report_metadata_html}')
os.makedirs(os.path.dirname(output_report_metadata_html), exist_ok=True)
with open(output_report_metadata_html, 'w', encoding='utf-8') as file:
    file.write("""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Report</title>
        <style>
            .table {
                border-collapse: collapse;
                width: 100%;
                font-family: 'Arial', sans-serif;
            }
            .table th,
            .table td {
                border: 1px solid #ddd;
                padding: 10px;
                text-align: left;
            }
            .table thead th {
                background-color: #f8f8f8;
                color: #333;
                font-weight: bold;
            }
            .table-striped tbody tr:nth-of-type(odd) {
                background-color: #f9f9f9;
            }
            .table-hover tbody tr:hover {
                background-color: #f1f1f1;
            }
        </style>
    </head>
    <body>
    """)
    file.write(html)
    file.write("""
    </body>
    </html>
    """)

###### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 