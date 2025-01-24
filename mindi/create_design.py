import argparse
from collections import defaultdict
from pathlib import Path
import pandas as pd
from termcolor import colored

def create_workflow_design(accession_list: str, tree: str) -> None:
    print(colored(f"Redirecting output to `design.csv`...", "green"))
    design = defaultdict(list)
    extract_id = lambda accession: '_'.join(Path(accession).name.split('_')[:2])
    with open(accession_list, mode="r", encoding="UTF-8") as f:
        for line in f:
            line = line.strip()
            accession_path = Path(line).resolve()
            accession_id = extract_id(accession_path)
            design["accession_id"].append(accession_id)
            design["accession"].append(accession_path)
    design = pd.DataFrame(design)
    design = merge_with_tree_of_life(design, tree)
    design.to_csv("design.csv", sep=",", mode="w", index=False, header=True)
    print(colored(f"Succesfully created `design.csv` file!", "green"))

def merge_with_tree_of_life(design, tree):
    tree_df = pd.read_csv(tree)[["#assembly_accession", 
                                   "organism_name", 
                                   "group", 
                                   "genome_size", 
                                   "gc_percent",
                                   "family",
                                   "order",
                                   "class",
                                   "phylum", 
                                   "kingdom", 
                                   "superkingdom"]]\
                    .set_index("#assembly_accession")
    merged_design = design.merge(tree_df, 
                                 how="left", 
                                 left_on="accession_id",
                                 right_on="#assembly_accession",
                                 right_index=True
                                 )
    # merged_design["genome_size"] = merged_design["genome_size"].astype(int)
    return merged_design
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Creates workflow design file.""")
    parser.add_argument("--accession", "-a", type=str, default="")
    parser.add_argument("--tree", "-t", type=str, default="tree_of_life_new.txt.gz")
    args = parser.parse_args()
    accession_list = args.accession
    tree = args.tree
    create_workflow_design(accession_list, tree=tree)
