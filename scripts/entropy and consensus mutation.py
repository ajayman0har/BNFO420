# -*- coding: utf-8 -*-
from Bio import SeqIO
import pandas as pd
from collections import Counter
import math

# Input files
aligned_fasta = "aligned.fasta"
wuhan_seq = ('MSDNGPQNQRNAPRITFGGPSDSTGSNQNGERSGARSKQRRPQGLPNNTASWFTALTQHGKEDLKFPRGQGVPINTNSSPDDQIGYYRRATRRIRGGDGKMKDLSPRWYFYYLGTGPEAGLPYGANKDGIIWVATEGALNTPKDHIGTRNPANNAAIVLQLPQGTTLPKGFYAEGSRGGSQASSRSSSRSRNSSRNSTPGSSRGTSPARMAGNGGDAALALLLLDRLNQLESKMSGKGQQQQGQTVTKKSAAEASKKPRQKRTATKAYNVTQAFGRRGPEQTQGNFGDQELIRQGTDYKHWPQIAQFAPSASAFFGMSRIGMEVTPSGTWLTYTGAIKLDDKDPNFKDQVILLNKHIDAYKTFPPTEPKKDKKKKADETQALPQRQKKQQTVTLLPAADLDDFSKQLQQSMSSADSTQA')



records = list(SeqIO.parse(aligned_fasta, "fasta")) #reading aligned sequences


seqs_22C = [str(rec.seq) for rec in records if "22C" in rec.description] #seperating out lineages
seqs_21J = [str(rec.seq) for rec in records if "21J" in rec.description]#using biopython to look at headers
all_seqs = [str(rec.seq) for rec in records]

def consensus(seqs): #create consensus sequences
    if not seqs:
        return "" #debugging
    consensus = []
    for i in range(len(seqs[0])):
        col = [s[i] for s in seqs if s[i] != "-"] #goes through each column for lineage
        consensus.append(Counter(col).most_common(1)[0][0] if col else "-") #counting the most common and just isolating the AA from the list of tuples
    return "".join(consensus)

def shannon_entropy(column): #shannon entropy method
    counts = Counter([res for res in column if res != "-"])
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((count / total) * math.log2(count / total) for count in counts.values())

def aa_frequencies(seqs):
    freq_table = []
    seq_len = len(seqs[0])
    for i in range(seq_len):
        col = [s[i] for s in seqs if s[i] != "-"]
        counts = Counter(col)
        total = sum(counts.values())
        freqs = {aa: round(count / total, 3) for aa, count in counts.items()} if total > 0 else {} #returning a list of dicts showing the amino acid frequencies at each position
        freq_table.append(freqs)
    return freq_table

# getting entropy per site
entropies = []
for i in range(len(all_seqs[0])):
    column = [seq[i] for seq in all_seqs]
    entropies.append(shannon_entropy(column))

#set consensus
cons_22C = consensus(seqs_22C)
cons_21J = consensus(seqs_21J)

# Compare to wuhan
freqs_22C = aa_frequencies(seqs_22C)
freqs_21J = aa_frequencies(seqs_21J)

# mutation threshold and comparison
freq_threshold = 0.05  # oonly in more than 5 percent of sequences
mutations = []

for i, ref_aa in enumerate(wuhan_seq, start=1):
    aa_22C = cons_22C[i - 1] if i - 1 < len(cons_22C) else "-"
    aa_21J = cons_21J[i - 1] if i - 1 < len(cons_21J) else "-"
    entropy = entropies[i - 1] if i - 1 < len(entropies) else 0.0

    freqs22 = freqs_22C[i - 1] if i - 1 < len(freqs_22C) else {}
    freqs21 = freqs_21J[i - 1] if i - 1 < len(freqs_21J) else {}

    # Include relevant amino acids are more than the frequency threshold
    relevant22 = {aa: f for aa, f in freqs22.items() if aa != ref_aa and f >= freq_threshold}
    relevant21 = {aa: f for aa, f in freqs21.items() if aa != ref_aa and f >= freq_threshold}

    if relevant22 or relevant21 or aa_22C != ref_aa or aa_21J != ref_aa:
        if aa_22C != ref_aa and aa_21J != ref_aa:
            mut_type = "Convergent"
            if aa_22C != aa_21J:
                mut_type += " (different AAs)"
        elif aa_22C != ref_aa:
            mut_type = "22C-specific"
        elif aa_21J != ref_aa:
            mut_type = "21J-specific"
        else:
            mut_type = "Minor-frequency mutation"

        mutations.append({"Position": i,"Wuhan_AA": ref_aa,"22C_AA": aa_22C, "21J_AA": aa_21J, "22C_Mutations": relevant22,"21J_Mutations": relevant21,"Type": mut_type, "Entropy": round(entropy, 4)})
#save results
df = pd.DataFrame(mutations)
df.to_csv("N_mutation_comparison_entropy.csv", index=False)

#full entropy of all locations
entropy_df = pd.DataFrame({"Position": list(range(1, len(entropies)+1)),"Entropy": entropies})
entropy_df.to_csv("N_entropy_profile.csv", index=False)

print("comparison complete")
df = pd.read_csv('N_mutation_comparison_entropy.csv')
df



