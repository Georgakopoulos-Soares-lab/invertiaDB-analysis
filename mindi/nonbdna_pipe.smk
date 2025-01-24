from termcolor import colored
import json
from collections import defaultdict
import csv
import numpy as np
import polars as pl
import pandas as pd
from mindi.minditool import MindiTool
from mindi.coverage.utils import parse_fasta, ProgressTracker
import threading
from mindi.scheduling import MiniBucketScheduler

buckets = config['buckets']
out = config['out']
pattern = config['pattern']
repr_pattern = '_'.join(c for c in pattern)
DESIGN = config['DESIGN']
merged_dest = Path(f"{out}/{repr_pattern}_merged").resolve()

def load_bucket(bucket: int) -> list[str]:
    global buckets
    global out
    print(f"Reading schedule data from '{out}/schedule_{buckets}.json'...")
    with open(f"{out}/schedule_{buckets}.json", mode="r", encoding="UTF-8") as f:
        return json.load(f)[str(bucket)]

def load_files(design: str) -> dict[str, str]:
    files = {}
    with open(design, mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=",")
        for row in reader:
            files[row['accession']] = row['accession_id']
    return files

rule all:
    input:
        expand("%s/%s_completed/%s_db_{mode}.parquet.snappy" % (out, repr_pattern, repr_pattern),
               mode=pattern),
        expand("%s/%s_completed/%s_db_{mode}.merged.parquet.snappy" % (out, repr_pattern, repr_pattern),
               mode=pattern)

rule schedule:
    input:
        DESIGN
    output:
        "%s/schedule_%s.json" % (out, config["buckets"]),
        "%s/nan_count.txt" % out,
        "%s/scheduled_%s.txt" % (out, repr_pattern),
    params:
        out=Path(config["out"]).resolve(),
        buckets=int(config["buckets"]),
        pattern=config['pattern'],
        suffix=config['suffix']
    run:
        assemblies = []
        with open(input[0], mode="r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter=",")
            for row in reader:
                assemblies.append(row['accession'])
        print(colored(f"Total assemblies detected: {len(assemblies)}.", "green"))
        if len(assemblies) == 0:
            raise ValueError(f'No assemblies were detected from the design file `{input[0]}`.')
        Path(params.out).mkdir(exist_ok=True)
        total_split = params.buckets
        mini_bucket_scheduler = MiniBucketScheduler()
        scheduled_files = mini_bucket_scheduler.schedule(assemblies, total_buckets=total_split)
        mini_bucket_scheduler.saveas(scheduled_files, output[0])
        
        # count nan values in each taxonomic rank
        ranks = ["family",
                 "order",
                 "class",
                 "phylum",
                 "kingdom",
                 "superkingdom"
                 ]
        design_df = pl.read_csv(input[0])
        total_nans = defaultdict(list)
        for rank in ranks:
            total_nans[rank] = design_df.select(
                                        pl.col(rank).is_null()
                                 )\
                                .sum()\
                                .item()

        with open(f"{out}/nan_count.txt", mode="w", encoding="utf-8") as f:
            for rank, nan_count in total_nans.items():
                print(nan_count)
                color = "red" if nan_count > 0 else "green"
                print(colored(f"NAN counts for taxonomic rank {rank}: {nan_count}.", color))
                f.write(f"{rank}\t{nan_count}\n")

        # create destination files; scheduled for delivery
        files = load_files(input[0])
        with open(f"{params.out}/scheduled_{repr_pattern}.txt", mode="w", encoding="utf-8") as f:
            for accession in assemblies:
               accession_name = Path(accession).name.split(f'.{params.suffix}')[0]
               f.write(files[accession] + "\t")
               for i, m in enumerate(pattern):
                   f.write(f"{params.out}/{repr_pattern}_extracted_accessions/{accession_name}_{m}.processed.tsv")
                   # append H-DNA extraction
                   if m == "MR":
                       f.write(f"\t{params.out}/{repr_pattern}_extracted_accessions/{accession_name}_H-DNA.processed.tsv")
                   if i < len(pattern) - 1:
                       f.write("\t")
                   else:
                       f.write("\n")

rule extractRepeats:
    input:
        DESIGN,
        '%s/schedule_%s.json' % (out, config["buckets"])
    output:
        touch("%s/%s_completed/bucket_{bucket}.%s.completed" % (out, repr_pattern, repr_pattern)),
        "%s/%s_completed/produced_merged_files_{bucket}.txt" % (out, repr_pattern),
    params:
        out=Path(config['out']).resolve(),
        pattern=config['pattern'],
        sleeping_time=float(config['log_sleep']),
    run:
        # PREPARE DESTINATION DIRECTORIES
        Path(f"{params.out}/{repr_pattern}_completed").mkdir(exist_ok=True, parents=True)
        merged_dest = Path(f"{out}/{repr_pattern}_merged").resolve()
        merged_dest.mkdir(exist_ok=True)

        # hash dir
        hash_dir = Path(f"{out}/.{repr_pattern}_hashes").resolve()
        hash_dir.mkdir(exist_ok=True)

        # temporary directory
        tempdir = Path(params.out).joinpath(f"{repr_pattern}_temp")
        print(f"Redirecting output to temporary directory --> `{tempdir}`.")

        # destination directory
        destination_dir = params.out.joinpath(f"{repr_pattern}_extracted_accessions")
        destination_dir.mkdir(exist_ok=True)

        # load files & prepare for extraction
        bucket = load_bucket(wildcards.bucket)
        mindi = MindiTool(tempdir=tempdir)
        total_accessions = len(bucket)
      
        # progress logging
        progress_tracking_log = params.out.joinpath("biologs", 
                                                    f"biolog_tracker_{repr_pattern}_{wildcards.bucket}.extraction.log")
        tracker = ProgressTracker(
                                 total_accessions=total_accessions,
                                 filename=progress_tracking_log,
                                 bucket_id=wildcards.bucket,
                                 sleeping_time=params.sleeping_time)
        logging_thread = threading.Thread(target=tracker.track_progress, daemon=True)
        logging_thread.start()

        # produced files
        produced_files = []
        produced_merged_files = []
        
        for accession in bucket:
            tracker.counter += 1
            print(f'Processing accession {accession}.')
            # destination = destination_dir.joinpath(accession_id + f".{mode}.csv")
            _ = mindi.extract(accession, 
                              pattern=params.pattern, 
                              destination_dir=destination_dir,
                              )
            accession_name = MindiTool.extract_name(accession)

            # validate extractions
            for m in params.pattern:
                # calculate md5hash
                extraction = mindi.fnp[m]
                os.system(f"md5sum {mindi.fnp[m]} >> {hash_dir}/md5_{wildcards.bucket}.{m}.hashes")

                if m == "STR" or m == "DR" or m == "IR" or m == "MR":
                    mindi.validate(accession, mode=m)
                    # remove file
                    os.remove(mindi.fn[m])
                
                # Extract H-DNA* from Mirror Repeats (mode=MR)
                if m == "MR":
                    df_hdna = mindi.extract_HDNA()
                    df_hdna.to_csv(
                               destination_dir.joinpath(accession_name + f"_H-DNA.processed.tsv"),
                               sep="\t",
                               mode="w",
                               index=False,
                               header=True
                               )
                    # Merge H-DNA extractions
                    if df_hdna.shape[0] > 0:
                        merged_bed_dest = merged_dest.joinpath(accession_name + f"_H-DNA.processed.merged.tsv")
                        df_hdna_merged = MindiTool.merge_frame(df_hdna[["seqID", "start", "end"]])\
                                                  .saveas(merged_bed_dest)
                        produced_merged_files.append(str(merged_bed_dest))

                # merge accession ? 
                merged_bed = mindi.merge(mode=m)
                if merged_bed is not None:
                    merged_bed_dest = merged_dest.joinpath(accession_name + f"_{m}.processed.merged.tsv")
                    merged_bed.saveas(merged_bed_dest)
                    produced_merged_files.append(str(merged_bed_dest))
            print(f'Accession {accession} has been processed succesfully.')
        with open(f"{out}/{repr_pattern}_completed/produced_merged_files_{wildcards.bucket}.txt", 
                  mode="w",
                  encoding="utf-8") as f:
            for merged_file in produced_merged_files:
                f.write(merged_file + "\n")
            
def load_delivery(params) -> dict[str, dict[str, str]]:
      delivered = {}
      patterns = params.pattern
      if "MR" in patterns:
          MR_idx = patterns.index("MR")
          patterns = patterns[:MR_idx+1] + ["H-DNA"] + patterns[MR_idx+1:]
      with open(f"{params.out}/scheduled_{repr_pattern}.txt", mode="r", encoding="utf-8") as f:
          data = csv.DictReader(f, 
                                delimiter="\t", 
                                fieldnames=["accession_id"] + [f"accession_name_{m}" for m in patterns]
                              )
          for row in data:
              delivered[row['accession_id']] = {m: row[f'accession_name_{m}'] for m in patterns}
      return delivered

def load_merged_delivery(params) -> dict[str, dict[str, str]]:
      delivered = load_delivery(params)
      delivered = {key: {m: merged_dest.joinpath(
                        Path(v).name\
                            .replace(".processed", ".processed.merged")
                            ) for m, v in value.items()} for key, value in delivered.items()}
      return delivered

rule reduceRepeats:
    input:
        DESIGN,
        expand("%s/%s_completed/bucket_{bucket}.%s.completed" % (out, repr_pattern, repr_pattern), 
               bucket=range(buckets))
    output:
        expand("%s/%s_completed/%s_db_{mode}.parquet.snappy" % (out, repr_pattern, repr_pattern),
               mode=pattern),
        expand("%s/%s_completed/%s_db_{mode}.merged.parquet.snappy" % (out, repr_pattern, repr_pattern),
               mode=pattern),
        "%s/%s_completed/empty_accessions.%s.txt" % (out, repr_pattern, repr_pattern),
        "%s/%s_completed/empty_accessions.enriched.%s.txt" % (out, repr_pattern, repr_pattern)
    params:
        out=Path(config["out"]).resolve(),
        buckets=config["buckets"],
        pattern=config['pattern'],
        step_threshold=float(config['step_threshold']),
    resources:
        mem_mb=60000
    run:
        target = params.out.joinpath(f"{repr_pattern}_extracted_accessions")
        print(f"Fetching extracted accessions from target {target}.")
        files = load_files(input[0])
        df_all = defaultdict(list)
        df_merged_all = defaultdict(list)
        empty_assemblies = defaultdict(list)

        # load expected extractions
        delivered = load_delivery(params)
        delivered_merged = load_merged_delivery(params)
        # load design
        design_df = pd.read_csv(input[0]) #, dtype={"genome_size": int})

        # hashdir
        hash_dir = Path(f"{out}/{repr_pattern}_hashes").resolve()
        hash_dir.mkdir(exist_ok=True)

        # patterns
        patterns = params.pattern
        if "MR" in patterns:
            MR_idx = patterns.index("MR")
            patterns = patterns[:MR_idx+1] + ["H-DNA"] + patterns[MR_idx+1:]

        # md5 hashes
        for m in patterns:
            os.system(f"cat {hash_dir}/md5_*.{m}.hashes > {hash_dir}/md5_all.{m}.hashes")
      
        # merge step / concatenating all extractions to a single dataframe
        # - for original
        # - and for merged (bed) files
        for file, file_id in files.items():
            print(f"Processing file `{file}`...")
            delivery = delivered[file_id]
            delivery_merged = delivered_merged[file_id]

            for m in patterns:
                # concatenate extractions
                assert m in delivery and m in delivery_merged, "Invalid parameter mapping."
                df = pd.read_table(delivery[m])
                if df.shape[0] == 0:
                    empty_assemblies[file].append(1)
                    print(f"mode {m} is EMPTY --> {delivery[m]}.")
                    continue
                empty_assemblies[file].append(0)
                df.loc[:, "#assembly_accession"] = file_id
                df_all[m].append(df)

                # concatenate merged extractions
                df_merged = pd.read_table(delivery_merged[m],
                                          header=None,
                                          names=["seqID", "start", "end", "counts"])
                df_merged.loc[:, "#assembly_accession"] = file_id
                df_merged_all[m].append(df_merged)

        # save empty accessions
        EMPTY = "%s/%s_completed/empty_accessions.%s.txt" % (out, repr_pattern, repr_pattern)
        EMPTY_ENRICHED = "%s/%s_completed/empty_accessions.enriched.%s.txt" % (out, repr_pattern, repr_pattern)
        with open(EMPTY, mode="w", encoding="UTF-8") as g:
            pattern_string = '\t'.join(m for m in patterns)
            g.write(f"#assembly_accession\taccession\t{pattern_string}\n")
            for file, values in empty_assemblies.items():
                g.write(files[file] + "\t" + file + "\t" + "\t".join(str(v) for v in values) + "\n")
        pd.read_table(EMPTY)\
                    .merge(design_df, 
                           right_on=["accession_id", "accession"],
                           left_on=["#assembly_accession", "accession"],
                           how="left"
                          )\
                    .drop(columns=['accession'])\
                    .to_csv(
                        EMPTY_ENRICHED,
                        sep="\t",
                        mode="w",
                        index=False,
                        header=True)

        # concatenate dataframes
        for mode in patterns:
            # concatenate original DB
            pd.concat(df_all[mode], axis=0)\
                .merge(design_df, 
                       left_on="#assembly_accession",
                       right_on="accession_id",
                       how="left"
                 )\
                .drop(columns=['accession'])\
                .to_parquet(f"{out}/{repr_pattern}_completed/{repr_pattern}_db_{mode}.parquet.snappy", 
                            compression="snappy",
                            engine="fastparquet")

            # concatenate & create merged DB
            pd.concat(df_merged_all[mode], axis=0)\
                .merge(design_df, 
                       left_on="#assembly_accession",
                       right_on="accession_id",
                       how="left"
                 )\
                .drop(columns=['accession'])\
                .to_parquet(f"{out}/{repr_pattern}_completed/{repr_pattern}_db_{mode}.merged.parquet.snappy", 
                            compression="snappy",
                            engine="fastparquet")

        # remove merged files for space efficiency
        print(colored(f"Removing intermediate merged bedfiles...", "blue"))
        for file, file_id in files.items():
            delivery_merged = delivered_merged[file_id]
            # remove merged accession for space efficiency
            for m in patterns:
                if Path(delivery_merged[m]).is_file():
                    pass
                    # os.remove(delivery_merged[m])
        print(colored(f"Intermediate merged bedfiles have been succesfully removed.", "green"))

        # calculate taxonomic densities
        for mode in patterns:
            print(colored(f"Calculating genomic densities for `{mode}`.", "blue"))
            if mode == "STR":
                seq_col = "sequence"
            else:
                seq_col = "sequenceOfArm"
            outfile = f"{out}/{repr_pattern}_completed/{repr_pattern}_db_{mode}.parquet.snappy"
            shell(f"""python calc_densities.py \
                        --dest {out} \
                        --source {outfile} \
                        --empty {EMPTY_ENRICHED} \
                        --mode {mode} \
                        --step_threshold {params.step_threshold} \
                        --seq_col {seq_col} \
                        --merged 0 \
                        --calculate_gc 1
                  """)
            print(colored(f"Genomic densities have succesfully been calculated for `{mode}`.", "green"))

            print(colored(f"Calculating merged genomic densities for `{mode}`.", "blue"))
            outfile = f"{out}/{repr_pattern}_completed/{repr_pattern}_db_{mode}.merged.parquet.snappy"
            shell(f"""python calc_densities.py \
                        --dest {out} \
                        --source {outfile} \
                        --empty {EMPTY_ENRICHED} \
                        --mode {mode} \
                        --step_threshold {params.step_threshold} \
                        --seq_col {seq_col} \
                        --merged 1 \
                        --calculate_gc 0
                  """)
            print(colored(f"Genomic densities for merged files have succesfully been calculated for `{mode}`.", "green"))

            # TODO
            # add biophysical properties calculation
        # pipeline exits
