team_member = {
    "name": "Sri Sathya Sandilya Garemilla",
    "slack_username": "@Garemilla Sandilya",
    "country": "United States",
    "hobby": "Coding",
    "affiliation": "Bioinformatics and Molecular Biochemistry",
    "favorite_gene": "FOS",
    "sequence_header": ">NC_000014.9:75278828-75282230 Homo sapiens chromosome 14, GRCh38.p14 Primary Assembly",
    "raw_sequence": """
AACCGCATCTGCAGCGAGCATCTGAGAAGCCAAGACTGAGCCGGCGGCCGCGGCGCAGCGAACGAGCAGT
GACCGTGCTCCTACCCAGCTCTGCTCCACAGCGCCCACCTGTCTCCGCCCCTCGGCCCCTCGCCCGGCTT
TGCCTAACCGCCACGATGATGTTCTCGGGCTTCAACGCAGACTACGAGGCGTCATCCTCCCGCTGCAGCA
GCGCGTCCCCGGCCGGGGATAGCCTCTCTTACTACCACTCACCCGCAGACTCCTTCTCCAGCATGGGCTC
GCCTGTCAACGCGCAGGTAAGGCTGGCTTCCCGTCGCCGCGGGGCCGGGGGCTTGGGGTCGCGGAGGAGG
AGACACCGGGCGGGACGCTCCAGTAGATGAGTAGGGGGCTCCCTTGTGCCTGGAGGGAGGCTGCCGTGGC
CGGAGCGGTGCCGGCTCGGGGGCTCGGGACTTGCTCTGAGCGCACGCACGCTTGCCATAGTAAGAATTGG
TTCCCCCTTCGGGAGGCAGGTTCGTTCTGAGCAACCTCTGGTCTGCACTCCAGGACGGATCTCTGACATT
AGCTGGAGCAGACGTGTCCCAAGCACAAACTCGCTAACTAGAGCCTGGCTTCTCCGGGGAGGTGGCAGAA
AGCGGCAATCCCCCCTCCCCCGGCAGCCTGGAGCACGGAGGAGGGATGAGGGAGGAGGGTGCAGCGGGCG
GGTGTGTAAGGCAGTTTCATTGATAAAAAGCGAGTTCATTCTGGAGACTCCGGAGCGGCGCCTGCGTCAG
CGCAGACGTCAGGGATATTTATAACAAACCCCCTTTCAAGCAAGTGATGCTGAAGGGATAACGGGAACGC
AGCGGCAGGATGGAAGAGACAGGCACTGCGCTGCGGAATGCCTGGGAGGAAAAGGGGGAGACCTTTCATC
CAGGATGAGGGACATTTAAGATGAAATGTCCGTGGCAGGATCGTTTCTCTTCACTGCTGCATGCGGCACT
GGGAACTCGCCCCACCTGTGTCCGGAACCTGCTCGCTCACGTCGGCTTTCCCCTTCTGTTTTGTTCTAGG
ACTTCTGCACGGACCTGGCCGTCTCCAGTGCCAACTTCATTCCCACGGTCACTGCCATCTCGACCAGTCC
GGACCTGCAGTGGCTGGTGCAGCCCGCCCTCGTCTCCTCCGTGGCCCCATCGCAGACCAGAGCCCCTCAC
CCTTTCGGAGTCCCCGCCCCCTCCGCTGGGGCTTACTCCAGGGCTGGCGTTGTGAAGACCATGACAGGAG
GCCGAGCGCAGAGCATTGGCAGGAGGGGCAAGGTGGAACAGGTGAGGAACTCTAGCGTACTCTTCCTGGG
AATGTGGGGGCTGGGTGGGAAGCAGCCCCGGAGATGCAGGAGCCCAGTACAGAGGATGAAGCCACTGATG
GGGCTGGCTGCACATCCGTAACTGGGAGCCCTGGCTCCAAGCCCATTCCATCCCAACTCAGACTCTGAGT
CTCACCCTAAGAAGTACTCTCATAGTTTCTTCCCTAAGTTTCTTACCGCATGCTTTCAGACTGGGCTCTT
CTTTGTTCTCTTGCTGAGGATCTTATTTTAAATGCAAGTCACACCTAGTCTGCAACTGCAGGTCAGAAAT
...
"""
}


def Sequence(header: str, raw_seq: str) -> str:
    lines = raw_seq.strip().splitlines()
    cleaned = "".join(line.strip() for line in lines if not line.startswith(">")).upper().replace("U", "T")
    return cleaned


def translate_dna(dna_seq: str, frame: int = 1, to_stop: bool = False) -> str:
    codon_table = {
        "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
        "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
        "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
        "TCT":"S","TCC":"S","TCA":"S","TCG":"S","AGT":"S","AGC":"S",
        "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
        "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
        "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
        "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
        "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
        "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
        "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
        "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
        "CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R",
        "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
    }
    if frame not in (1, 2, 3):
        raise ValueError("frame must be 1, 2, or 3")
    start = frame - 1
    protein = []
    for i in range(start, len(dna_seq)-2, 3):
        codon = dna_seq[i:i+3]
        aa = codon_table.get(codon, "X")
        if aa == "*" and to_stop:
            break
        protein.append(aa)
    return "".join(protein)



DNA = Sequence(team_member["sequence_header"], team_member["raw_sequence"])


protein_seq = translate_dna(DNA, frame=1)

# Print information
print("=== Team Member Information ===")
print(f"Name: {team_member['name']}")
print(f"Slack Username: {team_member['slack_username']}")
print(f"Country: {team_member['country']}")
print(f"Hobby: {team_member['hobby']}")
print(f"Affiliation: {team_member['affiliation']}")
print(f"Favorite Gene: {team_member['favorite_gene']}")
print(f"DNA Sequence Length: {len(DNA)} bp")

# Print protein
print(f"\n=== Protein Translation (Frame 1) ===\n{protein_seq[:300]}...") 
