import pandas as pd
import re

# File paths
ground_truth_file = 'ground_truth.csv'  # Adjust the path if necessary
input_file = 'input.csv'                # Adjust the path if necessary
output_file = 'output.csv'              # Adjust the output path if necessary

# Read input files with correct column names
ground_truth = pd.read_csv(ground_truth_file, sep=',', header=0, dtype=str)
input_data = pd.read_csv(input_file, sep=',', header=0, dtype=str)

# Function to parse gene_id into (Stx type, variant, allele, Accession)
def parse_gene_id(gene_id):
    if pd.isna(gene_id):
        return None, None, None, None
    match = re.match(r'(Stx\d+)([a-zA-Z]?)_([\d]+)_([\w.-]+)\|', gene_id)
    if match:
        return match.groups()  # (Stx type, Variant, Allele, Accession)
    return None, None, None, None

# Parse gene_id into structured columns for both input and ground truth
ground_truth[["stx_type", "stx_variant", "allele", "accession"]] = ground_truth["gene_id"].apply(lambda x: pd.Series(parse_gene_id(x)))
input_data[["stx_type", "stx_variant", "allele", "accession"]] = input_data["gene_id"].apply(lambda x: pd.Series(parse_gene_id(x)))

# Normalize accession numbers by replacing underscores with hyphens
ground_truth["accession"] = ground_truth["accession"].str.replace('_', '-')
input_data["accession"] = input_data["accession"].str.replace('_', '-')

# Prepare a dictionary to store ground truth data for each SeqID
ground_truth_dict = ground_truth.groupby("SeqID").agg(list).to_dict(orient='index')

# Initialize an empty list to store the output
output_rows = []

# Compare the input data with the ground truth data for each SeqID
for _, row in input_data.iterrows():
    seqid = row["SeqID"]
    input_stx_type, input_stx_variant, input_allele, input_accession = row["stx_type"], row["stx_variant"], row["allele"], row["accession"]

    if seqid in ground_truth_dict:
        truth_data = ground_truth_dict[seqid]
        truth_stx_type = truth_data["stx_type"]
        truth_stx_variant = truth_data["stx_variant"]
        truth_allele = truth_data["allele"]
        truth_accession = truth_data["accession"]
        strain = truth_data["Strain"][0]

        # Compare input vs truth
        stx_type_match = "T" if input_stx_type in truth_stx_type else "F"
        stx_variant_match = "T" if input_stx_variant in truth_stx_variant else "F"
        allele_match = "T" if input_allele in truth_allele else "F"
        accession_match = "T" if input_accession in truth_accession else "F"

        output_rows.append([strain, seqid, stx_type_match, stx_variant_match, allele_match, accession_match])
    else:
        # If no matching rows in the ground truth, mark as "F"
        output_rows.append(["-", seqid, "F", "F", "F", "F"])

# Convert results to DataFrame
output_df = pd.DataFrame(output_rows, columns=["Strain", "SeqID", "stx", "variant", "allele", "accession"])

# Compute percentage of fully matched rows (all True for stx_variant, allele, and accession)
accuracy = (output_df[["stx", "variant", "allele", "accession"]] == "T").all(axis=1).sum() / len(output_df) * 100

# Write output to file
output_df.to_csv(output_file, sep='\t', index=False)

# Append accuracy to the output file
with open(output_file, 'a') as f:
    f.write(f"\nAccuracy: {accuracy:.2f}%\n")

# Print accuracy
print(f"Accuracy: {accuracy:.2f}%")