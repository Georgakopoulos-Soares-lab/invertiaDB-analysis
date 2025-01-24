from collections import defaultdict, Counter

class DNA:

    def __init__(self, dna_string: str) -> None:
        if not isinstance(dna_string, str):
            raise TypeError(f"Invalid DNA string {dna_string}.")
        self.dna_string = dna_string

    def transmute(self, translation_table: int = 0, phase: int = 0) -> str:
        pass

    def __len__(self) -> int:
        return len(self.dna_string)
    
    @staticmethod
    def complement(nucleotide: str) -> str:
        if nucleotide == "a" or nucleotide == "A":
            return "t"
        if nucleotide == "t" or nucleotide == "T":
            return "a"
        if nucleotide == "c" or nucleotide == "C":
            return "g"
        if nucleotide == "g" or nucleotide == "G":
            return "c"
        if nucleotide == "n" or nucleotide == "N":
            return "n"
        raise ValueError(f"Invalid nucleotide {nucleotide}.")

    def count_kmers(self, kmer_length: int) -> dict[str, int]:
        if kmer_length == 1:
            return dict(Counter(self.dna_string))
        kmer_counts = defaultdict(list)
        for i in range(len(self) - kmer_length + 1):
            kmer_counts[self[i:i+kmer_length]] += 1
        return dict(kmer_counts)

    def __getitem__(self, key) -> str:
        if isinstance(key, slice):
            return self.dna_string[key.start: key.stop]
        if isinstance(key, int):
            return self.dna_string[key]
        raise ValueError()
    
    @staticmethod
    def reverse_complement(dna_string: str) -> str:
        return ''.join(DNA.complement(c) for c in dna_string)[::-1]

    def reverse_complement(self) -> str:
        return DNA.reverse_complement(self.dna_string)

    def _translation_table(self, translation_code: int) -> dict[str, int]:
        pass


class Protein:

    def __init__(self, peptide_string: str) -> None:
        if not isinstance(peptide_string, str):
            raise TypeError(f"Invalid Peptide string {peptide_string}.")
        self.peptide_string = peptide_string

    def biophysical(self):
        pass

    @classmethod
    def from_dna_string(cls, dna_string: str, translation_table: str, phase: int) -> "Protein":
        dna = DNA(dna_string)
        return cls(dna.transmute(translation_table, phase))


