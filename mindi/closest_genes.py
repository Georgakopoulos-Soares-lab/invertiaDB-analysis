import csv
from pybedtools import BedTool 
import argparse
import pandas as pd
from pathlib import Path

GFF_FIELDS = ["seqID", 
              "source", 
              "compartment", 
              "start", 
              "end", 
              "score", 
              "strand", 
              "phase", 
              "attributes"]
def read_gff(gff: str) -> pd.DataFrame:
    def _parse_attributes(metadata: str) -> dict[str, str]:
        attributes = {}
        for attribute in metadata.split(";"):
            key, value = attribute.split("=")
            attributes.update({key: value})
        return attributes
    def _get_attribute(attr: str, metadata: str) -> str:
        return _parse_attributes(metadata)[attr]
    df = pd.read_table(gff,
                         header=None,
                         names=GFF_FIELDS,
                         comment="#"
                         )
    df["start"] = df["start"] - 1
    df["gene_id"] = df["attributes"].apply(lambda metadata: _get_attribute("ID", metadata))
    df = df[df["compartment"].str.contains("gene")]
    df["biotype"] = df["attributes"].apply(lambda metadata: _get_attribute("gene_biotype", metadata))
    gff_bed = (
            BedTool.from_dataframe(df[["seqID", "start", "end", "gene_id", "biotype", "strand"]])
            .sort()
        )
    return gff_bed
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", type=str)
    parser.add_argument("--accession", type=str)
    args = parser.parse_args()
    gff = args.gff
    accession = args.accession
    gff_bed = read_gff(gff)
    extraction_bed = BedTool.from_dataframe(
            pd.read_table(accession, 
                          usecols=["chromosome", "start", "end"])
            ).sort()
    closest_genes = extraction_bed.closest(gff_bed, D='b')
    # closest_genes_df = pd.read_table(closest_genes.fn, header=None)
    closest_genes_df = []
    extract_id = lambda accession: "_".join(Path(accession).name.split("_")[:2])
    accession_id = extract_id(accession)
    with open(closest_genes.fn, 'r') as f:
        reader = csv.DictReader(f, 
                                delimiter="\t",
                                fieldnames=["seqID", 
                                            "start", 
                                            "end", 
                                            "chromosome", 
                                            "gene_start", 
                                            "gene_end", 
                                            "gene_id", 
                                            "gene_biotype", 
                                            "strand", 
                                            "distance"])
        for row in reader:
            start = int(row['start'])
            end = int(row['end'])
            gene_start = int(row['gene_start'])
            gene_end = int(row['gene_end'])
            gene_strand = row['strand']
            if gene_strand == '+':
                if end <= gene_start:
                    distance_from_TSS = gene_start - end + 1
                else:
                    distance_from_TSS = max(0, start - gene_start)
                if end <= gene_end:
                    distance_from_TES = gene_end - end 
                else:
                    distance_from_TES = max(0, start - gene_end + 1)
            else:
                if end <= gene_start:
                    distance_from_TES = gene_start - end + 1
                else:
                    distance_from_TES = max(0, start - gene_start)
                if end <= gene_end:
                    distance_from_TSS = gene_end - end 
                else:
                    distance_from_TSS = max(0, start - gene_end + 1)
            closest_genes_df.append(
                    row | {"distance_from_TSS": distance_from_TSS,
                           "distance_from_TES": distance_from_TES})
    closest_genes_df = pd.DataFrame(closest_genes_df)
    closest_genes_df.to_csv(f"{accession_id}_g4.closest_gene.txt", 
                            sep="\t", 
                            mode="w")
    breakpoint()
