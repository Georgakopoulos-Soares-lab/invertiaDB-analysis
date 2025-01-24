# Mindi: A Non B-DNA extractor tool for data analysis

__author__ = "Nikol Chantzi"
__version__ = "1.0.1"
__email__ = "nmc6088@psu.edu"

import sys
import os
import shutil
import gzip
import logging
import uuid
import csv
import re
from termcolor import colored
from pathlib import Path
import json
import tempfile
import subprocess
from dotenv import load_dotenv
from typing import ClassVar, Optional, Iterator
import pandas as pd
import polars as pl
from pybedtools import BedTool
from mindi.coverage.utils import parse_fasta
# from mindi.tailhunter import hunt_tail

# config
csv.field_size_limit(sys.maxsize)
logging.basicConfig(level=logging.WARNING, 
                    format="%(asctime)s:%(levelname)s:%(message)s")

class MindiTool:

    REGEX_FIELDS: ClassVar[list[str]] = [
                    "seqID",
                    "start",
                    "end",
                    "sequence",
                    "strand",
                    "gc_content",
                    "length",
                    "stackerLength",
                    "spacerLength",
                    ]
    TAIL_FIELDS: ClassVar[list[str]] = ["seqID", 
                                        "start", 
                                        "end", 
                                        "sequence", 
                                        "length", 
                                        "consensus", 
                                        "sru", 
                                        "consensus_repeats"]
    HDNA_max_at_content: ClassVar[float] = 0.8
    HDNA_min_pyrine: ClassVar[float] = 0.9
    HDNA_max_pyrine: ClassVar[float] = 0.9
    # please do not modify; makes program insanely slow
    max_DR_rep: ClassVar[int] = 2_000

    def __init__(self, tempdir: Optional[str] = None,
                       nonBDNA: Optional[str] = None) -> None:
        load_dotenv()
        if nonBDNA is None:
            nonBDNA = os.getenv('nonBDNA')
        if tempdir is None:
            self.tempdir = Path().cwd()
        else:
            self.tempdir = Path(tempdir).resolve()
            self.tempdir.mkdir(exist_ok=True)
        self.nonBDNA = nonBDNA
        self.fn: dict = {}
        self.fnp: dict = {}
        
        with open(Path(__file__).resolve().parent.joinpath("defaults.json"),
                  mode="r",
                  encoding="utf-8") as f:
            data = json.load(f)
        self.defaults = data["params"]
        self.MINDI_FIELDS = data["headers"]
        if not Path(self.nonBDNA).is_file():
            raise FileNotFoundError(f"Invalid executable path {self.nonBDNA}.")

    @staticmethod
    def extract_name(accession: str) -> str:
        accession = Path(accession).resolve()
        accession_name = accession.name
        if accession_name.endswith('.fna.gz'):
            return accession_name.split('.fna.gz')[0]
        if accession_name.endswith('.fa.gz'):
            return accession_name.split('.fa.gz')[0]
        if accession_name.endswith('.fna'):
            return accession_name.split('.fna')[0]
        if accession_name.endswith('.fa'):
            return accession_name.split('.fa')[0]
        return accession_name
        
    @staticmethod
    def extract_id(accession: str) -> str:
        return '_'.join(Path(accession).name.split("_")[:2])

    def _generate_repeats(self, accession: str, pattern: list[str], **kwargs) -> "MindiTool":
        self.reset()
        accession = Path(accession).resolve()
        if not accession.is_file():
            raise FileNotFoundError(f'Unable to detect accession {accession}.')
        accession_name = MindiTool.extract_name(accession)
        # accession_tmp_dir = self.tempdir.joinpath(accession_name)
        accession_tmp_dir = tempfile.TemporaryDirectory(dir=self.tempdir,
                                                        prefix=f"{accession_name}_{'_'.join(pattern)}.",
                                                        )
        accession_tmp_dir_path = self.tempdir.joinpath(accession_tmp_dir.name)
        # if accession_tmp_dir.is_dir():
        #    shutil.rmtree(accession_tmp_dir)
        cur_dir = os.getcwd()
        # handle zipped files
        if accession.name.endswith(".gz"):
            # gzip compression strategy
            with tempfile.NamedTemporaryFile(dir=accession_tmp_dir_path,
                                             delete=False,
                                             suffix='.fna') as unzipped_tmp:
                with gzip.open(accession, 'rb') as handler:
                    unzipped_tmp.write(handler.read())
                accession = Path(unzipped_tmp.name).name
        os.chdir(accession_tmp_dir_path)
        rand_accession_name = str(uuid.uuid4())
        # unfortunately non b-gfa breaks with large names 
        if len(accession.name) > 60:
            shutil.copy(accession, accession.name)
            accession = accession.name
        command = f"{self.nonBDNA} -seq {accession} -out {rand_accession_name} -skipAPR -skipSlipped -skipCruciform -skipTriplex -skipWGET"
        seen_modes = set()
        skipped = set()
        for m in pattern:
            if m == "IR":
                command += f" -minIRrep {kwargs['IR']['minrep']} -maxIRspacer {kwargs['IR']['maxspacer']}"
            elif m == "MR":
                command += f" -minMRrep {kwargs['MR']['minrep']} -maxMRspacer {kwargs['MR']['maxspacer']}"
            elif m == "DR":
                command += f" -minDRrep {kwargs['DR']['minrep']} -maxDRrep {kwargs['DR']['maxrep']} -maxDRspacer {kwargs['DR']['maxspacer']}"
            seen_modes.add(m)
        for m in ["IR", "MR", "DR", "Z", "GQ", "STR"]:
            if m not in seen_modes:
                if m == "IR":
                    command += " -skipIR"
                elif m == "MR":
                    command += " -skipMR"
                elif m == "DR":
                    command += " -skipDR"
                elif m == "STR":
                    command += " -skipSTR"
                elif m == "Z":
                    command += " -skipZ"
                elif m == "GQ":
                    command += " -skipGQ"
                skipped.add(m)
        if not seen_modes:
            raise ValueError(f'Unknown pattern match `{pattern}`.')
        print(f"Input pattern detected: `{pattern}`.\nRunning command `{command}`...")
        _ = subprocess.run(command, shell=True,
                                    check=True,
                                    stdout=subprocess.DEVNULL,
                                    # stderr=subprocess.DEVNULL,
                            )
        # check if operation was succesful
        for mode in pattern:
            if not Path(rand_accession_name + f"_{mode}.tsv").is_file():
                raise FileNotFoundError(f"Failed to extract {mode} for {accession}.")
        # continue with processing
        for mode in pattern:
            shutil.move(rand_accession_name + f"_{mode}.tsv", accession_name + f"_{mode}.tsv")
            shutil.move(rand_accession_name + f"_{mode}.gff", accession_name + f"_{mode}.gff")
            destination = self.tempdir.joinpath(accession_name + f'_{mode}.tsv')
            if destination.is_file():
                os.remove(destination)
            out_tsv = accession_tmp_dir_path.joinpath(accession_name + f'_{mode}.tsv')
            out_gff = accession_tmp_dir_path.joinpath(accession_name + f'_{mode}.gff')
            shutil.move(out_tsv, destination)
            # remove redundant files
            os.unlink(out_gff)
            # modify tmp file pointer
            self.fn[mode] = self.tempdir.joinpath(accession_name + f'_{mode}.tsv')
        accession_tmp_dir.cleanup()
        os.chdir(cur_dir)
        return self

    @property
    def extractions(self) -> list[str]:
        return list(self.fn.keys())

    def reset(self) -> None:
        self.fn = {}
        self.fnp = {}

    def moveto(self, dest: str) -> None:
        dest = Path(dest).resolve()
        for m in self.extractions:
            shutil.move(self.fnp[m], dest)
            self.fnp[m] = dest.joinpath(self.fnp[m].name)
            if m == "IR" or m == "DR" or m == "MR":
                shutil.move(self.fnp[m], dest)
            else:
                shutil.move(self.fn[m], dest)

    def to_dataframe(self, mode: str, usecols: bool = True) -> pd.DataFrame:
        # STILL EXPERIMENTAL
        if mode == 'RE' or mode == "TAIL":
            tmp_file = self.fn[mode]
        else:
            tmp_file = self.fnp[mode]
        if usecols:
            mindi_frame = pd.read_table(tmp_file)
            if mode == "IR" or mode == "DR" or mode == "MR":
                mindi_frame.loc[:, "sequenceOfSpacer"] = mindi_frame["sequenceOfSpacer"].fillna(".")
        else:
            mindi_frame = pd.read_table(tmp_file, header=None, skiprows=1)
        return mindi_frame

    def set_tempdir(self, tempdir: os.PathLike[str]) -> None:
        print(f"Resetting tempdir to --> `{tempdir}`.")
        self.tempdir = Path(tempdir).resolve()
        self.tempdir.mkdir(exist_ok=True)

    def extract(self, accession: str, pattern: list[str] | str, destination_dir: str, **kwargs):
        # filter duplicates
        if isinstance(pattern, str):
            pattern = list(set(pattern.split("|")))
        if not isinstance(pattern, list):
            raise TypeError(f"Unknown type {type(pattern)}.")

        # process arguments/parameters for each mode
        for m in pattern:
            if m == 'IR' or m == 'DR' or m == 'MR':
                if m in kwargs:
                    kwargs[m]['minrep'] = kwargs[m].get('minrep', self.defaults[m]['minrep'])
                    kwargs[m]['maxspacer'] = kwargs[m].get('maxspacer', self.defaults[m]['maxspacer'])
                else:
                    kwargs[m] = {}
                    for attr, val in self.defaults[m].items():
                        kwargs[m]['maxspacer'] = self.defaults[m]['maxspacer']
                        kwargs[m]['minrep'] = self.defaults[m]['minrep']
            if m == 'DR':
                if m in kwargs:
                    kwargs[m]['maxrep'] = kwargs[m].get('maxrep', self.defaults[m]['maxrep'])
                else:
                    kwargs[m]['maxrep'] = self.defaults[m]['maxrep']
            if m == 'Z' or m == 'STR' or m == 'GQ':
                kwargs[m] = {}

        labels = {"MR": "Mirror Repeats",
                  "IR": "Inverted Repeats",
                  "DR": "Direct Repeats",
                  "STR": "Short Tandem Repeats",
                  "GQ": "GQuadruplex",
                  "APR": "Poly-A Tract",
                  "RE": "G4 (Consensus Motif)",
                  "Z": "Z-DNA",
                  }
        stars = "=" * 10
        for m in pattern:
           print(f"{stars}<{labels[m]}> parameters{stars}")
           for attr, val in kwargs[m].items():
               print(f"{attr}={val}")
           print(stars)

        self._generate_repeats(accession=accession, pattern=pattern, **kwargs)
        for mode in pattern:
            if mode == "IR" or mode == "DR" or mode == "MR" or mode == "STR":
                self.process_table(mode=mode, dest_path=destination_dir, **kwargs[mode])
        return self
    
    @staticmethod
    def get_HDNA(mindi_table: pd.DataFrame) -> pd.DataFrame:
        mindi_table = mindi_table.copy()
        mindi_table.loc[:, "pyrine"] = (mindi_table["sequence"].str.count("a|g")).div(mindi_table["sequenceLength"])
        mindi_table.loc[:, "pyrimidine"] = 1.0 - mindi_table["pyrine"]
        # mindi_table.loc[:, "pyrimidine"] = (mindi_table["sequence"].str.count("c|t")).div(mindi_table["sequenceLength"])
        mindi_table.loc[:, "at_content"] = (mindi_table["sequence"].str.count("a|t")).div(mindi_table["sequenceLength"])
        mindi_table = mindi_table[(mindi_table["at_content"] <= 0.8) & ((mindi_table["pyrimidine"] >= 0.9) | (mindi_table["pyrine"] >= 0.9))]\
                            .reset_index(drop=True)
        mindi_table.loc[:, "motif_strand"] = (mindi_table["pyrine"] >= 0.5).astype(int).apply(lambda x: "+" if x == 1 else "-")
        # H-DNA Template
        # non template --> ag rich strand
        # template --> ct rich strand
        return mindi_table

    def extract_HDNA(self) -> pd.DataFrame:
        # assert current extraction is of mode 'Mirror Repeat'
        if 'MR' not in self.fnp:
            raise ValueError("Cannot proceed with H-DNA filtering. Mirror Repeats have not been extracted.")
        mindi_table = pd.read_table(self.fnp["MR"])
        return MindiTool.get_HDNA(mindi_table)

    def extract_g4(self, accession: os.PathLike[str],
                            stacker: str = "g",
                            minrep: int = 3,
                            multiplicity: int = 3) -> "MindiTool":
        regex = RegexExtractor(
                               stacker=stacker,
                               minrep=minrep,
                               multiplicity=multiplicity
                               )
        accession_name = MindiTool.extract_name(accession)
        accession_tmp_dir = self.tempdir.joinpath(accession_name)
        accession_tmp_dir.mkdir(exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=accession_tmp_dir,
                                         prefix=accession_name + "_RE.",
                                         delete=False, 
                                         suffix=".tsv", 
                                         mode="w") as file:
            dict_writer = csv.DictWriter(file, 
                                         delimiter="\t", 
                                         fieldnames=MindiTool.REGEX_FIELDS)
            dict_writer.writeheader()
            for table in regex.parse_g4(accession):
                for row in table:
                    dict_writer.writerow(row)
            self.fn["RE"] = file.name
        return self

    def cleanup(self) -> None:
        if self.fnp and Path(self.fnp).is_file():
            os.remove(self.fnp)
        else:
            print(colored(f"WARNING! Processed file {self.fnp} does not exist.", "red"))

    def extract_tails(self, accession: os.PathLike[str], minrepeat: int = 8) -> "MindiTool":
        with tempfile.NamedTemporaryFile(dir=self.tempdir,
                                         delete=False,
                                         suffix=".tail",
                                         mode="w") as file:
            dict_writer = csv.DictWriter(
                                         file,
                                         delimiter="\t",
                                         fieldnames=MindiTool.TAIL_FIELDS
                                         )
            dict_writer.writeheader()
            for seqID, sequence in parse_fasta(accession):
                # TODO implmenet hunt tail < import  
                for mononucleotide_tail in hunt_tail(sequence=sequence,
                                                     minrepeat=minrepeat,
                                                     seqID=seqID
                                                     ):
                    dict_writer.writerow(mononucleotide_tail)
            self.fn["TAIL"] = file.name
        return self

    @staticmethod
    def complement(nucleotide: str) -> str:
        if nucleotide == "a":
            return "t"
        if nucleotide == "t":
            return "a"
        if nucleotide == "g":
            return "c"
        if nucleotide == "c":
            return "g"
        if nucleotide == "n":
            return "n"
        raise ValueError(f"Unknown nucleotide {nucleotide}.")

    @staticmethod
    def reverse(kmer: str) -> str:
        return ''.join(MindiTool.complement(c) for c in kmer)[::-1]

    def validate(self, accession: os.PathLike[str], mode: str) -> None:
        if mode == "Z" or mode == "GQ":
            raise NotImplementedError(f"Validation for extraction type {mode} has yet to be implemented.")
        print(colored(f"Validating accession `{accession}` for {mode=}...", "blue"))
        df = self.to_dataframe(mode=mode)
        for seqID, seq in parse_fasta(accession):
            temp = df[df['seqID'] == seqID]
            total = temp.shape[0]
            validated = 0
            for _, row in temp.iterrows():
                start = row['start']
                end = row['end']
                sequence = row['sequence']
                arm_seq = row['sequenceOfArm']
                spacer_seq = row['sequenceOfSpacer']
                if spacer_seq == ".":
                    spacer_seq = ""
                fasta_seq = seq[start: end]
                # general assertions
                assert sequence == fasta_seq
                if mode == "MR":
                    assert len(arm_seq) * 2 + len(spacer_seq) == len(sequence) == end - start
                    assert arm_seq + spacer_seq + arm_seq[::-1] == fasta_seq
                elif mode == "DR":
                    assert len(arm_seq) * 2 + len(spacer_seq) == len(sequence) == end - start
                    assert arm_seq + spacer_seq + arm_seq == fasta_seq 
                elif mode == "IR":
                    assert len(arm_seq) * 2 + len(spacer_seq) == len(sequence) == end - start
                    assert arm_seq + spacer_seq + MindiTool.reverse(arm_seq) == fasta_seq
                elif mode == "STR":
                    consensus_repeats = row['consensus_repeats']
                    sru = row['sru']
                    assert arm_seq * consensus_repeats == fasta_seq 
                    assert sru * consensus_repeats == len(fasta_seq)
                validated += 1
            assert validated == total
        print("OK!")
        print(colored(f'Accession {accession} has passed all checks.', 'green'))

    @staticmethod
    def merge_frame(df: pd.DataFrame) -> BedTool:
        return BedTool.from_dataframe(df)\
                    .sort()\
                    .merge(c="3", o="count")

    def merge(self, mode: str) -> BedTool:
        df_extraction = pd.read_table(self.fnp[mode], 
                                      usecols=["seqID", "start", "end"])
        if df_extraction.shape[0] == 0:
            return
        return BedTool.from_dataframe(df_extraction)\
                        .sort()\
                        .merge(c="3", o="count")

    def merged_df(self, mode: str) -> pd.DataFrame:
        bed_merged = self.merge(mode=mode)
        if bed_merged is None:
            return
        return pd.read_table(
                        bed_merged.fn,
                        header=None,
                        names=["seqID", "start", "end", "overlapping_count"]
                    )\
                .assign(total_bp=lambda ds: ds['end'] - ds['start'])

    def read(self, mode: str) -> pl.DataFrame:
        return pl.read_csv(self.fn[mode], separator="\t")

    def readp(self, mode: str) -> pl.DataFrame:
        return pl.read_csv(self.fnp[mode], separator="\t")

    def process_table(self, mode: str,
                            minrep: Optional[int] = None,
                            maxrep: Optional[int] = None,
                            maxspacer: Optional[int] = None,
                            dest_path: Optional[str] = None
                            ) -> "MindiTool":
        to_drop = ["Source",
                   "Type",
                   "Score",
                   "Strand",
                   "Subset",
                   "Permutations",
                   "Sequence",
                   "Start",
                   "Stop"]
        if dest_path:
            dest_path = Path(dest_path).resolve()
        tmp_file = tempfile.NamedTemporaryFile(dir=self.tempdir,
                                               prefix=str(self.fn[mode]).replace(".tsv", "") + ".",
                                               suffix=".processed.tsv",
                                               delete=False,
                                               mode="w"
                                            )
        nucleotides = {'a', 'g', 'c', 't'}
        FIELDS = self.MINDI_FIELDS.copy()
        # add STR specific fields
        if mode == "STR":
            FIELDS += ["Repeated", "sru", "consensus_repeats"]

        with tmp_file as fout:
            fout_writer = csv.DictWriter(fout, delimiter="\t", fieldnames=FIELDS)
            fout_writer.writeheader()
            with open(self.fn[mode], mode="r", encoding="UTF-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    arm_length = int(row['Repeat'])
                    spacer_length = int(row['Spacer'])
                    sequence = row['Sequence']
                    # check validity of the entry
                    if isinstance(minrep, int) and arm_length < minrep:
                        logging.warning(f"Invalid record detected with arm length {arm_length} lesser than the minimum arm length of {minrep}.")
                        continue
                    elif isinstance(maxspacer, int) and spacer_length > maxspacer:
                        logging.warning(f"Invalid record detected with spacer length {maxspacer} greater than the maximum spacer length of {maxspacer}.")
                        continue
                    elif isinstance(maxrep, int) and arm_length > maxrep:
                        logging.warning(f"Invalid record detected with arm length {maxrep} greater than the maximum arm length of {maxrep}.")
                        continue

                    start = int(row['Start']) - 1
                    end = int(row['Stop'])
                    sequence_length = int(row['Length'])
                    total_coordinate_length = end - start
                    if sequence_length < total_coordinate_length:
                        # sequence = sequence[:sequence_length]
                        # end = end - (total_coordinate_length - sequence_length)
                        print(colored(f'Invalid sequence length detected for {self.fn} on (start,end)=({start},{end}) with sequence {sequence}.', 'red'))
                        raise ValueError(f'Invalid sequence length detected for {self.fn}.')
                    # find sequence of arm
                    repeat = int(row['Repeat'])
                    sequence_of_arm = sequence[:repeat]
                    # find spacer
                    del row['Spacer']
                    if mode == "IR" or mode == "DR" or mode == "MR":
                        true_spacer_length = sequence_length - 2 * repeat
                        right_arm = sequence[repeat+true_spacer_length:]
                        if any(n not in nucleotides for n in right_arm) or any(n not in nucleotides for n in sequence_of_arm):
                            logging.warning("Invalid record detected with non-AGCT nucleotide present in sequence of arm.")
                            continue
                        # spacer = sequence[repeat:repeat+spacer_length]
                        true_spacer = sequence[repeat:repeat+true_spacer_length]
                        if len(true_spacer) == 0:
                            true_spacer = "."
                        # skip maximum spacer length
                        if (isinstance(maxspacer, int) and true_spacer_length > maxspacer):
                            # invalid record?
                            raise ValueError(f"Invalid record detected. Max spacer was found {true_spacer_length} but cannot exceed {maxspacer}.")
                    else:
                        # set null spacer for all other NON BDNA modes
                        true_spacer = "."
                        true_spacer_length = 0

                    # remove excess in STR repeats
                    if mode == "STR":
                        sru = len(sequence_of_arm)
                        consensus_repeats, remainder = row['Repeated'][1:].split('+')
                        consensus_repeats = int(consensus_repeats)
                        remainder = int(remainder)
                        if remainder == sru:
                            sequence = sequence[:sequence_length-sru]
                            remainder = 0
                        elif sequence_length%sru == 0:
                            consensus_repeats += remainder // sru
                            remainder = 0
                        assert remainder == sequence_length%sru
                        sequence = sequence[:sequence_length - sequence_length%sru]
                        sequence_length = len(sequence)
                        assert consensus_repeats * sru  == sequence_length
                        row["sru"] = sru 
                        row["consensus_repeats"] = consensus_repeats
                        # change end 
                        end = start + sequence_length

                    # process composition
                    composition = re.search(r"(\d+)A/(\d+)C/(\d+)G/(\d+)T", row['Composition'])
                    a_content = composition.group(1)
                    c_content = composition.group(2)
                    g_content = composition.group(3)
                    t_content = composition.group(4)
                    row.update({
                            "start": start,
                            "end": end,
                            "sequenceOfArm": sequence_of_arm,
                            "sequenceOfSpacer": true_spacer,
                            "spacer": true_spacer_length,
                            "sequence": sequence,
                            "arm_a": a_content,
                            "arm_c": c_content,
                            "arm_g": g_content,
                            "arm_t": t_content,
                        })
                    for col in to_drop:
                        if col in row:
                            del row[col]
                    row["sequenceLength"] = row.pop("Length")
                    row["seqID"] = row.pop("Sequence_name")
                    row["armLength"] = row.pop("Repeat")
                    row["spacerLength"] = row.pop("spacer")
                    row["composition"] = row.pop("Composition")
                    fout_writer.writerow(row)

        if dest_path:
            target_name = dest_path.joinpath(self.fn[mode].name.replace(".tsv", ".processed.tsv"))
            shutil.move(tmp_file.name, target_name)
            self.fnp[mode] = target_name
        else:
            self.fnp[mode] = tmp_file.name
        return self

# Regex Pattern Extraction Implementations
# Does not report overlapping patterns 

class RegexExtractor:

    def __init__(self, stacker: str = "g", multiplicity: int = 3, minrep: int = 3) -> None:
        self.stacker = stacker
        self.multiplicity = multiplicity
        self.complementary_stacker = RegexExtractor.reverse(self.stacker)
        self.minrep = minrep
        self.motifs = [
                       "%s{%s,}[agct]{1,7}" % (self.stacker, self.minrep) * self.multiplicity + "%s{%s,}" % (self.stacker, self.minrep),
                       "%s{%s,}[agct]{1,7}" % (self.complementary_stacker, self.minrep) * self.multiplicity + "%s{%s,}" % ( self.complementary_stacker, self.minrep)
                       ]

        self.strict_motif = "%s" % (self.stacker * self.minrep) + "[agct]{1,7}"
        self.strict_motif = self.strict_motif * self.multiplicity + self.strict_motif
        # self.complementary_strict_motif

    @staticmethod
    def complement(nucleotide: str) -> str:
        # please do not use heap allocation for no reason =)
        if nucleotide == "a":
            return "t" 
        if nucleotide == "t":
            return "a"
        if nucleotide == "g":
            return "c"
        if nucleotide == "c":
            return "g"
        if nucleotide == "n":
            return "n"
        raise ValueError(f'Invalid nucleotide {nucleotide}.')

    @staticmethod
    def reverse(kmer: str) -> str:
        return ''.join(RegexExtractor.complement(n) for n in kmer)[::-1]

    @staticmethod
    def _extract_tandem(seq: str, minrep: int) -> Iterator[dict]:
        pattern = re.compile(r"([AGCT]{%s,})\1+" % minrep)
        matches = re.finditer(pattern, seq)
        for match in matches:
            start = match.start()
            end = match.end()
            consensus_motif = match.group(1)
            sru = len(consensus_motif)
            sequence = match.group(0)
            sequence_length = len(sequence)
            consensus_repeats = sequence_length // sru
            yield {
                    "start": start,
                    "end": end,
                    "sru": sru,
                    "consensus_motif": consensus_motif,
                    "consensus_repeats": consensus_repeats,
                    "sequence": sequence,
                    "length": sequence_length,
                    }

    @staticmethod
    def extract_tandem(seq: str, minrep: int, seqID: Optional[str] = None) -> pd.DataFrame:
        tandem_df = []
        for row in RegexExtractor._extract_tandem(seq=seq, minrep=minrep):
            tandem_df.append(row)
        tandem_df = pd.DataFrame(tandem_df)
        if seqID:
            tandem_df["seqID"] = seqID
        return tandem_df
    
    @staticmethod
    def _extract_direct(seq: str, minrep: int, minspacer: int, maxspacer: int) -> Iterator[dict]:
        pattern = re.compile(r"([AGCT]{%s,})([AGCT]{%s,%s}?)\1" % (minrep, minspacer, maxspacer))
        matches = re.finditer(pattern, seq)
        for match in matches:
            start = match.start()
            end = match.end()
            sequence_of_arm = match.group(1)
            sequence_of_spacer = match.group(2)
            arm_length = len(sequence_of_arm)
            spacer_length = len(sequence_of_spacer)
            sequence = match.group(0)
            sequence_length = len(sequence)
            yield {
                    "start": start,
                    "end": end,
                    "sequence": sequence,
                    "length": sequence_length,
                    "sequence_of_arm": sequence_of_arm,
                    "sequence_of_spacer": sequence_of_spacer,
                    "arm_length": arm_length,
                    "spacer_length": spacer_length,
                    }

    @staticmethod
    def extract_direct(seq: str, minrep: int, minspacer: int, maxspacer: int, seqID: Optional[str] = None) -> pd.DataFrame:
        direct_df = []
        for row in RegexExtractor._extract_direct(seq=seq, 
                                                  minrep=minrep, 
                                                  minspacer=minspacer, 
                                                  maxspacer=maxspacer):
            direct_df.append(row)
        direct_df = pd.DataFrame(direct_df)
        if seqID:
            direct_df["seqID"] = seqID
        return direct_df

    def parse_g4(self, accession: str | os.PathLike[str]) -> Iterator[list[dict]]:
        for seqID, seq in parse_fasta(accession):
            for motif in self.motifs:
                matches = re.finditer(motif, seq)
                table = []
                for match in matches:
                    sequence = match.group()
                    gc_content = sequence.count('g') + sequence.count('c')
                    strand = "+" if sequence.startswith(self.stacker) else "-"
                    if strand == "+":
                        stacker = self.stacker
                    else:
                        stacker = self.complementary_stacker
                    stackers = re.findall("%s{%s,}" % (stacker, self.minrep), sequence)
                    total_stacker_len = sum(map(len, stackers))
                    spacers = re.sub("%s{%s,}" % (stacker, self.minrep), "|", sequence)
                    total_spacer_len = sum(map(len, spacers.split("|")))
                    start = match.start()
                    end = match.end()
                    assert seq[start: end] == sequence, f"Invalid sequence detected at chromosome {seqID} at position ({start},{end})."
                    table.append({
                            "seqID": seqID,
                            "start": start,
                            "end": end,
                            "sequence": sequence,
                            "strand": strand,
                            "gc_content": gc_content,
                            "length": len(sequence),
                            "stackerLength": total_stacker_len,
                            "spacerLength": total_spacer_len,
                        })
                table = sorted(
                               table,
                               key=lambda x: (x['seqID'], x['start']),
                               reverse=False
                               )
                yield table

    def extract_g4(self, accession: str) -> pd.DataFrame:
        g4_df = []
        for table in self.parse_g4(accession=accession):
            g4_df.extend(table)
        g4_df = pd.DataFrame(g4_df)
        return g4_df
