import math
import itertools
import os
import sys

# STEP 1 : Generates  all possible 'm+1 mers'
def get_m_plus_1_mers(m): 
    bases = ['A', 'C', 'G', 'T']
    return ["".join(p) for p in itertools.product(bases, repeat=m+1)]

# STEP 2 : Train a Markov model on the full set of sequences provided
def build_model(sequences, m):
    n = m + 1
    #Initialize counts with pseudocount = 1 
    m_plus_1_mers = get_m_plus_1_mers(m)
    counts = {nmer: 1 for nmer in m_plus_1_mers}
    
    #Fill counts using sliding window logic 
    for seq in sequences:
        for i in range(len(seq) - n + 1):
            nmer = seq[i : i + n]
            if nmer in counts:
                counts[nmer] += 1
    
    # STEP 3: Building the Transition Probability Matrix 
    prefixes = sorted(list(set(nmer[:m] for nmer in m_plus_1_mers)))
    transitions = {}
    bases_list = ['A', 'C', 'G', 'T']
    
    for prefix in prefixes:
        # Sum of all possible transitions from this prefix 
        row_sum = sum(counts[prefix + b] for b in bases_list)
        # Calculate P(base | prefix) 
        transitions[prefix] = {b: counts[prefix + b] / row_sum for b in bases_list}
        
    return transitions

# STEP 4 : Taking inputs from the user
fasta_path = input("Enter the name of your FASTA file (Ex : abc.fasta). Ensure that it is in the same directory as the code : ")
m_input = input("Enter the order of the Markov model (m): ")

try:
    m = int(m_input)
except ValueError:
    print("Error: Markov order must be an integer.")
    sys.exit()

if not os.path.exists(fasta_path):
    print(f"Error: Could not find {fasta_path}. Ensure it is in the directory.")
    sys.exit()

sequences = []
with open(fasta_path, "r") as f:
    for line in f:
        clean_line = line.strip().upper()
        # Skip empty lines or header lines starting with '>'
        if clean_line and not clean_line.startswith(">"):
            sequences.append(clean_line)

if not sequences:
    print("Error: No DNA sequences found in the file.")
    sys.exit()

n = m + 1

# STEP 5: Training using all sequences in the file
print(f"Training Markov model of order {m} on {len(sequences)} sequences...")
model_matrix = build_model(sequences, m)

# STEP 6: Scoring each sequence and saving the output in a text file
output_filename = f"log_likelihood_scores_for_order{m}.txt"
print(f"Calculating scores for each line and saving to {output_filename}...")

with open(output_filename, "w") as f:
    for seq in sequences:
        log_p = 0
        # If the sequence is shorter than n, it cannot be scored by an m-order model
        if len(seq) >= n: 
            for i in range(len(seq) - n + 1):
                nmer = seq[i : i + n]
                prefix, suffix = nmer[:m], nmer[m:]
                prob = model_matrix[prefix][suffix]
                log_p += math.log(prob)
            
            f.write(f"{log_p:.6f}\n")
        else:
            # Sequences too short for the chosen order m
            f.write("0.000 (Sequence too short)\n")

print(f"Successfully processed {len(sequences)} lines. Results saved.")
