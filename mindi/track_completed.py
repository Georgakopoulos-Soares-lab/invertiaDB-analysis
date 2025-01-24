from mindi.coverage.utils import load_bucket 
import pandas as pd
from termcolor import colored
from pathlib import Path
import json

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""""")
    parser.add_argument("--dir", type=str, default="")
    parser.add_argument("--suffix", type=str, default="processed.tsv")
    parser.add_argument("--schedule", type=str, default="repeat_out/schedule_2.json")
    parser.add_argument("--design", type=str, default="design.csv")

    args = parser.parse_args()
    suffix = args.suffix
    design = args.design
    schedule = args.schedule
    dir = Path(args.dir).resolve()

    # read schedule 
    with open(schedule, mode="r") as f:
        buckets = json.load(f)
    # read design
    # design = set(pd.read_csv(design, usecols=["accession_id"]))
        
    total_buckets = len(buckets)
    extract_id = lambda accession: '_'.join(Path(accession).name.split('_')[:2]) if Path(accession).name.count("_") > 1 else Path(accession).name.split(".")[0]
    extracted_files = {extract_id(file) for file in dir.glob(f"*.{suffix}")}
    
    expected_buckets = {bucket_id: set([extract_id(file) for file in bucket]) for bucket_id, bucket in buckets.items()}
    retrieved_buckets = {bucket_id: set() for bucket_id in buckets}
    for bucket_id, bucket in buckets.items():
        for file in bucket:
            accession_id = extract_id(file)
            if accession_id in extracted_files:
                retrieved_buckets[bucket_id].add(accession_id)

    for bucket_id in buckets:
        expected_buckets[bucket_id] = expected_buckets[bucket_id] - retrieved_buckets[bucket_id]
    
    for bucket_id, bucket in expected_buckets.items():
        if len(bucket) == 0:
            print(colored(f"Bucket {bucket_id} has been extracted.", "green"))
        else:
            print(colored(f"Bucket {bucket_id} is not extracted.", "red"))
            print("Remaining files:")
            for file in bucket:
                print(file)
