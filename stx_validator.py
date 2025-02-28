import pandas as pd
import re

# File paths
ground_truth_file = 'ground_truth.csv'  # Adjust the path if necessary
input_file = 'kma.csv'                # Adjust the path if necessary
output_file = 'output.csv'                  # Adjust the output path if necessary

# Read input files with correct column names
ground_truth = pd.read_csv(ground_truth_file, sep=',', header=0, dtype=str)
input_data = pd.read_csv(input_file, sep=',', header=0, dtype=str)

# Function to parse gene_id into (Stx type, variant, allele, Accession)
def parse_gene_id(gene_id):
    if pd.isna(gene_id):
        return None, None, None, None
    match = re.match(r'(Stx\d+)([a-zA-Z]?)_([\d]+)_([\w.-]+)', gene_id)
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

# Initialize an empty list to store the output and sets to track negative controls and missing SeqIDs
output_rows = []
negative_controls = set()
not_found_seqids = set()
missing_in_input = set(ground_truth["SeqID"].unique())
false_positives = []
new_seqids = set()
processed_seqids = set()
stx_mismatch_seqids = set()

# Compare the input data with the ground truth data for each SeqID
for _, row in input_data.iterrows():
    seqid = row["SeqID"]
    input_gene_id = row["gene_id"]

    if seqid in ground_truth_dict:
        truth_data = ground_truth_dict[seqid]
        strain = truth_data["Strain"][0]

        # Check for negative control
        if "Neg_control_no_vtx" in truth_data["gene_id"]:
            negative_controls.add(seqid)
            continue

        # Remove found SeqID from missing_in_input set
        missing_in_input.discard(seqid)

        # Initialize stx_type_match
        stx_type_match = "F"

        # Find the matching gene_id from ground truth
        ground_truth_geneid = next((gid for gid in truth_data["gene_id"] if gid.replace('_', '-') == input_gene_id.replace('_', '-')), None)

        if ground_truth_geneid:
            # Compare input vs truth
            stx_type_match = "T" if row["stx_type"] in truth_data["stx_type"] else "F"
            stx_variant_match = "T" if row["stx_variant"] in truth_data["stx_variant"] else "F"
            allele_match = "T" if row["allele"] in truth_data["allele"] else "F"
            accession_match = "T" if row["accession"] in truth_data["accession"] else "F"

            output_rows.append([strain, seqid, ground_truth_geneid, input_gene_id, stx_type_match, stx_variant_match, allele_match, accession_match])

            # Track SeqIDs where stx_type is "F"
            if stx_type_match == "F":
                stx_mismatch_seqids.add(seqid)
        else:
            # If no exact match, find the closest match
            closest_match = None
            for gid in truth_data["gene_id"]:
                if gid.split('_')[0] == input_gene_id.split('_')[0]:  # Match on Stx type
                    closest_match = gid
                    break
            if closest_match:
                stx_type_match = "T" if row["stx_type"] in truth_data["stx_type"] else "F"
                stx_variant_match = "T" if row["stx_variant"] in truth_data["stx_variant"] else "F"
                allele_match = "T" if row["allele"] in truth_data["allele"] else "F"
                accession_match = "F"
                output_rows.append([strain, seqid, closest_match, input_gene_id, stx_type_match, stx_variant_match, allele_match, accession_match])
                # Track SeqIDs where stx_type is "F"
                if stx_type_match == "F":
                    stx_mismatch_seqids.add(seqid)
            else:
                stx_mismatch_seqids.add(seqid)  # Add to stx_mismatch_seqids if no closest match is found
    else:
        # Track SeqIDs not found in the ground truth
        not_found_seqids.add(seqid)

# Ensure only one match is printed for SeqIDs present only once in the ground truth
unique_seqids = ground_truth["SeqID"].value_counts()
for seqid, count in unique_seqids.items():
    if count == 1:
        output_rows = [row for i, row in enumerate(output_rows) if row[1] != seqid or i == next(i for i, r in enumerate(output_rows) if r[1] == seqid)]

# Convert output list to DataFrame
output_df = pd.DataFrame(output_rows, columns=["Strain", "SeqID", "ground_truth_geneid", "input_gene_id", "stx", "variant", "allele", "accession"])

# Ensure unique entries for each SeqID and input_gene_id
output_df = output_df.drop_duplicates(subset=["SeqID", "input_gene_id", "ground_truth_geneid", "stx", "variant", "allele", "accession"])

# Compute percentage of fully matched rows (all True for stx, variant, allele, and accession)
overall_accuracy = (output_df[["stx", "variant", "allele", "accession"]] == "T").all(axis=1).mean() * 100

# Compute accuracy for each column
stx_accuracy = (output_df["stx"] == "T").mean() * 100
variant_accuracy = (output_df["variant"] == "T").mean() * 100
allele_accuracy = (output_df["allele"] == "T").mean() * 100
accession_accuracy = (output_df["accession"] == "T").mean() * 100

# Exclude negative controls from missing_in_input
missing_in_input -= negative_controls
# Exclude SeqIDs containing "Neg_control_no_vtx" from missing_in_input
missing_in_input = {seqid for seqid in missing_in_input if not any("Neg_control_no_vtx" in gene_id for gene_id in ground_truth_dict[seqid]["gene_id"])}

# Write cleaned output to CSV
output_df.to_csv(output_file, sep='\t', index=False, columns=["Strain", "SeqID", "ground_truth_geneid", "input_gene_id", "stx", "variant", "allele", "accession"])

# Append accuracy and negative control information to the output file
with open(output_file, 'a') as f:
    f.write(f"\nOverall Accuracy: {overall_accuracy:.2f}%\n")
    f.write(f"stx Accuracy: {stx_accuracy:.2f}%\n")
    f.write(f"variant Accuracy: {variant_accuracy:.2f}%\n")
    f.write(f"allele Accuracy: {allele_accuracy:.2f}%\n")
    f.write(f"accession Accuracy: {accession_accuracy:.2f}%\n")
    if negative_controls:
        f.write("\nSeqID(s) identified as negative controls (not vtx):\n")
        for seqid in negative_controls:
            f.write(f"{seqid}\n")
    if missing_in_input:
        f.write("\nSeqID(s) from ground truth file are absent in the user input file:\n")
        for seqid in missing_in_input:
            f.write(f"{seqid}\n")
    if not_found_seqids:
        f.write("\nSeqID(s) from input file are not found in the ground truth data (false positives):\n")
        for seqid in not_found_seqids:
            f.write(f"{seqid}\n")
    if stx_mismatch_seqids:
        for seqid in stx_mismatch_seqids:
            if seqid not in output_df["SeqID"].values:
                f.write(f"\nstx_type is not matching and its a false one didn't use it in our accuracy calculation:{seqid}\n")

# Print accuracy
print(f"Overall Accuracy: {overall_accuracy:.2f}%")
print(f"stx Accuracy: {stx_accuracy:.2f}%")
print(f"variant Accuracy: {variant_accuracy:.2f}%")
print(f"allele Accuracy: {allele_accuracy:.2f}%")
print(f"accession Accuracy: {accession_accuracy:.2f}%")
if negative_controls:
    print("\nSeqID(s) identified as negative controls (not vtx):")
    for seqid in negative_controls:
        print(seqid)
# Print SeqIDs missing in the input file
if missing_in_input:
    print("\nSeqID(s) from ground truth absent in the user input file:")
    for seqid in missing_in_input:
        print(seqid)
# Print SeqIDs not found in the ground truth
if not_found_seqids:
    print("\nSeqID(s) not found in the ground truth data (false positives):")
    for seqid in not_found_seqids:
        print(seqid)
# Print SeqIDs where stx_type is "F" and no match was found
if stx_mismatch_seqids:
    for seqid in stx_mismatch_seqids:
        if seqid not in output_df["SeqID"].values:
            print("\nstx_type is not matching and its a false one didn't use it in our accuracy calculation:" + seqid)
