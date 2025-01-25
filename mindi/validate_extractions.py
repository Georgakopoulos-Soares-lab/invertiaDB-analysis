from utils import parse_fasta, load_bucket
import argparse
from pathlib import Path
import pandas as pd

def find_right_arm(sequence_of_arm, mode):
    if mode == "IR":
        return "".join(
                        {
                            "a": "t", 
                            "t": "a",
                            "g": "c",
                            "c": "g"
                            }[n] for n in sequence_of_arm)[::-1]
    elif mode == "MR":
        return sequence_of_arm[::-1]
    elif mode == "DR":
        return sequence_of_arm
    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""""")
    parser.add_argument("--bucket_id", type=int, default=0)
    parser.add_argument("--schedule_path", type=str, default="schedule_100.json")
    parser.add_argument("--extractions_path", type=str, default="repeat_out")
    parser.add_argument("--mode", type=str, choices=["IR", "MR", "DR"], default="IR")

    args = parser.parse_args()
    bucket_id = args.bucket_id
    schedule = args.schedule_path
    extractions_path = Path(args.extractions_path).resolve()
    mode = args.mode

    files = load_bucket(bucket_id, schedule)
    
    extract_id = lambda accession: '_'.join(Path(accession).name.split("_")[:2])
    extract_id = lambda accession: Path(accession).name.split('.')[0]
    extractions = {extract_id(accession): accession for accession in extractions_path.glob("*.csv")}
    
    total_accessions = len(files)
    print(f"Initializing validation procedure for bucket {bucket_id} (mode {mode}). Total: {total_accessions} accessions.")
    
    total_ok = 0
    flag = 0
    total_non_empty = 0
    for file in files:
        # print(f"Processing accession {file}...")

        accession_id = extract_id(file)
        extracted_accession = extractions.get(accession_id)

        if extracted_accession is None:
            print(f"Skipping validation for file {file} since no extracted file was detected.")
            continue

        
        extracted_df = pd.read_table(extracted_accession)


        for seqID, seq in parse_fasta(file):
            
            temp = extracted_df[extracted_df['seqID'] == seqID]
            total_found = temp.shape[0]
            total_non_empty += int(total_found > 0)
            total_validated = 0

            for _, row in temp.iterrows():
                flag = 1
                start = int(row['start'])
                end = int(row['end'])
                sequence_of_arm = row['sequenceOfArm']
                sequence_of_spacer = row['sequenceOfSpacer']
                sequence = row['sequence']
                sequence_length = int(row['sequenceLength'])
                arm_length = int(row['armLength'])
                spacer_length = int(row['spacerLength'])

                if sequence_of_spacer == ".":
                    sequence_of_spacer = ""

                section = seq[start: end]
                
                right_arm = find_right_arm(sequence_of_arm, mode)
                assert section == sequence == sequence_of_arm + sequence_of_spacer + right_arm
                assert len(section) == sequence_length == arm_length * 2 + spacer_length
                total_validated += 1

            assert total_validated == total_found
        
        total_ok += 1

        # print(f"Accession {file} has passed all validations!")
    
    assert flag == 1
    print(f"Bucket: {bucket_id}. Total {total_ok} accessions have been validated out of {total_accessions}! Total non empty {total_non_empty}.")

                


            






