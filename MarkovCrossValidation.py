# This is a Markov model based classifier that predicts whether a given transcription factor binds to a given DNA sequence.
# Input : 1. m = order of the markov model (takes values from 0 to 10)
#         2. k = number of folds in the cross-validation (takes values from 3 to 5)
#         3. tsv_file = tsv file for a particular chromosome which contains start and stop positions of the 200bp bins
#         4. TF = name of the transcription factor of interest
# Output : 1. k ROC curves
#          2. k precision-recall curves
#          3. average area under the k ROC and precision-recall curves

import pandas as pd
from Bio import SeqIO
import os
import random
import itertools
import math
import time
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import matplotlib.pyplot as plt
import sys

# --- START TOTAL TIMER ---
start_total_time = time.time()

# Step 1 : Obtain the FASTA sequences of the bins in the tsv file

tsv_path = input("Enter the name of your TSV file (e.g., chr1_200bp_bins.tsv). Ensure it is in the working directory : ")

chrom_num = tsv_path.split('_')[0]
fasta_path = f"{chrom_num}.fasta"

if not os.path.exists(fasta_path):
    print(f"Error: Could not find {fasta_path}. Please ensure it is in the working directory.")
    sys.exit()
else:

    print(f"Loading {chrom_num}.fasta")
    record = SeqIO.read(fasta_path, "fasta")
    
    df = pd.read_csv(tsv_path, sep='\t')
    
    bins = {}
    for index, row in df.iterrows():
        start = int(row['start'])
        end = int(row['end'])
        seq_segment = str(record.seq[start:end])
        bin_id = f"bin{index + 1}"
        bins[bin_id] = seq_segment

    print(f"Successfully extracted {len(bins)} bins.")

# Step 2 : Seggregation of the bins into 'bound' and 'unbound' groups based on the TF of interest

tf_options = ['CTCF', 'REST', 'EP300']

print("Choice of TFs : CTCF, REST and EP300")
tf_name = input("Enter TF name : ").strip().upper()

if tf_name not in tf_options:
    print("Error: Invalid TF. Choose CTCF, REST, or EP300.")
    sys.exit()
else:
    def separate_bins(df, tf_col, bins):
        bound_seqs = []
        unbound_seqs = []
        
        for index, row in df.iterrows():
            bin_id = f"bin{index + 1}"
            seq = bins[bin_id].upper()
            
            if row[tf_col] == 'B':
                bound_seqs.append(seq)
            elif row[tf_col] == 'U':
                unbound_seqs.append(seq)
        
        print(f"Number of Bound bins: {len(bound_seqs)}")
        print(f"Number of Unbound bins: {len(unbound_seqs)}")
        return bound_seqs, unbound_seqs
    
    bound_seqs, unbound_seqs = separate_bins(df, tf_name, bins)

# Step 3 : Randomly dividing the bins into 'k' buckets such that each bucket gets roughly equal number of bins

k_folds = int(input("Enter the number of folds for cross-validation (must lie between 3-5): "))

random.shuffle(bound_seqs) 
random.shuffle(unbound_seqs)

def split_into_buckets(lst, k):

    folds = []
    for i in range(k):
        folds.append([]) 
    
    for index in range(len(lst)):
        bucket_index = index % k
        folds[bucket_index].append(lst[index])
        
    return folds

bound_folds = split_into_buckets(bound_seqs, k_folds)
unbound_folds = split_into_buckets(unbound_seqs, k_folds)

print(f"Shuffled and split {len(bound_seqs)} bound bins into {k_folds} folds.")
print(f"Shuffled and split {len(unbound_seqs)} unbound bins into {k_folds} folds.")

# Step 4 : List which contains all possible 'm+1 mers'.

m_input = (input("Enter the order of the Markov model (m) between 0 to 10 : "))
m = int(m_input)
n = m + 1

bases = ['A', 'T', 'G', 'C']

m_plus_1_mers = []

for p in itertools.product(bases, repeat=n):
    word = "".join(p)
    m_plus_1_mers.append(word)

print(f"Total number of m+1 mers generated : {len(m_plus_1_mers)}")

# Step 5 : Counting all the occurances of m+1 mers in the bins

def fill_counts(sequences, target_dict, n):
    for seq in sequences:
        for i in range(len(seq) - n + 1):
            current_nmer = seq[i : i + n]
            if current_nmer in target_dict:
                target_dict[current_nmer] = target_dict[current_nmer] + 1

# Step 6 : Building the Transition Probability Matrix

def transition_prob_matrix(counts_dict, m):

    prefixes = []
    for nmer in counts_dict.keys():    
        p = nmer[:m]
        if p not in prefixes:
            prefixes.append(p)
    prefixes.sort()

    transitions = {}
    bases = ['A', 'C', 'G', 'T']
    
    for prefix in prefixes:
        row_sum = 0 
        for b in bases:
            row_sum += counts_dict[prefix + b]  
        probabilities = {}
        for b in bases:
            probabilities[b] = counts_dict[prefix + b] / row_sum    
            
        transitions[prefix] = probabilities
        
    return transitions

# Step 7 : Training the model on k-1 buckets and using the remaining 1 for cross-validation

all_fold_test_labels = []
all_fold_test_scores = []

for i in range(k_folds):
    print(f"Training and Scoring Fold {i+1}...")
    
    training_data_bound = []
    training_data_unbound = []
    
    for j in range(k_folds):
        if j != i:      
            training_data_bound.extend(bound_folds[j])      
            training_data_unbound.extend(unbound_folds[j])
            
    b_counts = {}
    u_counts = {}
    
    for nmer in m_plus_1_mers:
        b_counts[nmer] = 1     
        u_counts[nmer] = 1
        
    fill_counts(training_data_bound, b_counts, n)   
    fill_counts(training_data_unbound, u_counts, n)
    
    b_matrix = transition_prob_matrix(b_counts, m)  
    u_matrix = transition_prob_matrix(u_counts, m)
    
    test_sequences = bound_folds[i] + unbound_folds[i]      
    test_labels = [1] * len(bound_folds[i]) + [0] * len(unbound_folds[i])   
    
    fold_scores = []
    for seq in test_sequences:
        log_p_b = 0
        log_p_u = 0

        for idx in range(len(seq) - n + 1):
            nmer = seq[idx : idx + n]
            prefix = nmer[:m]
            suffix = nmer[m:]
            
            log_p_b = log_p_b + math.log(b_matrix[prefix][suffix])
            log_p_u = log_p_u + math.log(u_matrix[prefix][suffix])
            
        fold_scores.append(log_p_b - log_p_u)   
    
    all_fold_test_labels.append(test_labels)
    all_fold_test_scores.append(fold_scores)

print("Cross-validation scoring is complete")

# Step 8 : Constructing the ROC and Precision-Recall Curves for all 'k' folds

roc_title = f"ROC Curves for {k_folds} fold validation for order {m}"
pr_title = f"Precision-Recall Curves for {k_folds} fold validation for order {m}"

# --- PLOT 1: ROC CURVES ---

plt.figure(figsize=(8, 6), dpi=150)
roc_aucs = []

for i in range(k_folds):
    fpr, tpr, _ = roc_curve(all_fold_test_labels[i], all_fold_test_scores[i])
    fold_auc = auc(fpr, tpr)
    roc_aucs.append(fold_auc)
    plt.plot(fpr, tpr, label=f'Fold {i+1} (AUC = {round(fold_auc, 3)})')

plt.plot([0, 1], [0, 1], 'k--', label='Random Guess')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title(roc_title)
plt.legend(loc='lower right')
plt.grid(alpha=0.3)
roc_filename = f"ROC_Curves_for_{k_folds}_fold_validation_for_order_{m}.png"
plt.savefig(roc_filename, bbox_inches='tight', dpi=150)
plt.show()


# --- PLOT 2: PRECISION-RECALL CURVES ---

plt.figure(figsize=(8, 6), dpi=150)
pr_aucs = []

for i in range(k_folds):
    precision, recall, _ = precision_recall_curve(all_fold_test_labels[i], all_fold_test_scores[i])
    fold_pr_auc = auc(recall, precision)
    pr_aucs.append(fold_pr_auc)
    plt.plot(recall, precision, label=f'Fold {i+1} (AUC = {round(fold_pr_auc, 3)})')

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title(pr_title)
plt.legend(loc='upper right')
plt.grid(alpha=0.3)
pr_filename = f"Precision-Recall_Curves_for_{k_folds}_fold_validation_for_order_{m}.png"
plt.savefig(pr_filename, bbox_inches='tight', dpi=150)
plt.show()

# --- END TOTAL TIMER AND LOG ---
end_total_time = time.time()
total_execution_duration = end_total_time - start_total_time

print("-" * 30)
print(f"TIME LOG FOR {k_folds} FOLDS for MM ORDER {m}:")
print(f"Total Execution Time: {total_execution_duration:.4f} seconds")
print("-" * 30)
