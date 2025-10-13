from Bio import Entrez, SeqIO
import pandas as pd
import time
import random
Entrez.email = "manohara@vcu.edu"
search_term = (
    "Severe acute respiratory syndrome coronavirus 2[Organism] "
    "AND human[host] AND 29000:31000[Sequence Length]"
) #searches for sars cov 2 sequqneces that are between 29k and 31k bp

handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1000) #queries ncbi and will return the ids of the sequences
record = Entrez.read(handle)
handle.close()

id_list = record["IdList"] # this will hold all the NCBI accession numbers
print(f"Found {len(id_list)} genomes — sampling 10 for N protein extraction.\n")


sample_ids = random.sample(id_list, 10) #randomly sample 10 from the id list


n_proteins = [] #fasta format
metadata = [] # id, oprganism, date , protein length

for seq_id in sample_ids:
    try:
        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text") #downloading the full genbank record and retype = gb will include the annotations
        record = SeqIO.read(handle, "genbank") # parses the record into a python object
        handle.close()

        for feature in record.features:
            if feature.type == "CDS":
                # Check for Nucleocapsid protein
                product = feature.qualifiers.get("product", [""])[0].lower()
                if "nucleocapsid" in product:
                    protein_seq = feature.qualifiers.get("translation", [""])[0]
                    n_proteins.append(f">{seq_id} {record.description}\n{protein_seq}\n")

                    metadata.append({
                        "id": seq_id,
                        "organism": record.annotations.get("organism", ""),
                        "collection_date": record.annotations.get("date", ""),
                        "protein_length": len(protein_seq)
                    })
                    break  # Stop after finding the N protein
        time.sleep(0.3) #this is to avoid overloading the server, 0.3 seconds between each request
    except Exception as e:
        print(f"Skipped {seq_id}: {e}") 


with open("sarscov2_N_proteins.fasta", "w") as f:
    for entry in n_proteins:
        f.write(entry)

df = pd.DataFrame(metadata)
df.to_csv("sarscov2_N_metadata.csv", index=False)
print(" saved: sarscov2_N_proteins.fasta")
print(" saved: sarscov2_N_metadata.csv")

