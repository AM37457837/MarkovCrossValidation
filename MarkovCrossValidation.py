# This is a Markov model based classifier that predicts whether a given transcription factor binds to a given DNA sequence.
# Input : 1. m = order of the markov model (takes values from 0 to 10)
#         2. k = number of folds in the cross-validation (takes values from 3 to 5)
#         3. tsv_file = tsv file for a particular chromosome which contains start and stop positions of the 200bp bins
#         4. TF = name of the transcription factor of interest
# Output : 1. k ROC curves
#          2. k precision-recall curves
#          3. average area under the k ROC and precision-recall curves

# Step 1 : Obtain the FASTA sequences of the bins in the tsv file
# Please install 'biopython' and 'pandas' libraries

import pandas as pd
from Bio import SeqIO
import os

# Input = TSV file
tsv_path = input("Enter the name of your TSV file (e.g., chr1_200bp_bins.tsv). Ensure it is in the working directory : ")

# Make sure to name the fasta file which contains the DNA sequence of the chromosome as {chrom_num}.fasta (Ex : chr15.fasta)
chrom_num = tsv_path.split('_')[0]
fasta_path = f"{chrom_num}.fasta"

if not os.path.exists(fasta_path):
    print(f"Error: Could not find {fasta_path}. Please ensure it is in the working directory.")
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

import random

k_folds = int(input("Enter the number of folds for cross-validation (must lie between 3-5): "))

random.shuffle(bound_seqs) # randomly shuffling the bins
random.shuffle(unbound_seqs)

def split_into_buckets(lst, k):

    folds = []
    for i in range(k):
        folds.append([]) # creating k empty buckets
    
    for index in range(len(lst)):
        bucket_index = index % k
        folds[bucket_index].append(lst[index])
        
    return folds

bound_folds = split_into_buckets(bound_seqs, k_folds)
unbound_folds = split_into_buckets(unbound_seqs, k_folds)

print(f"Shuffled and split {len(bound_seqs)} bound bins into {k_folds} buckets.")
print(f"Shuffled and split {len(unbound_seqs)} unbound bins into {k_folds} buckets.")

# Step 3 : List which contains all possible 'm+1 mers'.
# Since we are interested in a 'mth' order markov model, we want the probabability of getting a nucleotide, given 'm' previous nucleotides. Hence we need all possible 'm+1' mers

import itertools

m_input = (input("Enter the order of the Markov model (m) between 0 to 10 : "))
m = int(m_input)
n = m + 1

bases = ['A', 'T', 'G', 'C']

m_plus_1_mers = []

for p in itertools.product(bases, repeat=n):
    word = "".join(p)
    m_plus_1_mers.append(word)


print(f"Total number of m+1 mers generated : {len(m_plus_1_mers)}")


# Step 4 : Counting all the occurances of m+1 mers in the bins

def fill_counts(sequences, target_dict, n):
    for seq in sequences:
        for i in range(len(seq) - n + 1):
            current_nmer = seq[i : i + n]
            if current_nmer in target_dict:
                target_dict[current_nmer] = target_dict[current_nmer] + 1

# Step 5 : Building the Transition Probability Matrix

def transition_prob_matrix(counts_dict, m):

    prefixes = []
    for nmer in counts_dict.keys():    # getting the first 'm' charachters (the prefix) of every m+1 mer
        p = nmer[:m]
        if p not in prefixes:
            prefixes.append(p)
    prefixes.sort()

    transitions = {}
    bases = ['A', 'C', 'G', 'T']
    
    for prefix in prefixes:
        row_sum = 0 
        for b in bases:
            row_sum += counts_dict[prefix + b]  # getting all possible m+1 mers from the prefix, which is an m mer
        probabilities = {}
        for b in bases:
            probabilities[b] = counts_dict[prefix + b] / row_sum    # Since we added a pseudocount of 1, this sum will always be at least 4, preventing "division by zero" error.
            
        transitions[prefix] = probabilities
        
    return transitions

    # Example : for m=2, Transition Probability of AA --> G = P('AAG') = P(G|AA) = freq (AAG) / Sum of freq (AAA, AAC, AAG, AAT)
                                                
# Step 6 : Training the model on k-1 buckets and using the remaining 1 for cross-validation

import math

all_fold_test_labels = []
all_fold_test_scores = []

for i in range(k_folds):
    print(f"Training and Scoring Fold {i+1}...")
    
    # 1. Training Data (k-1 folds)
    training_data_bound = []
    training_data_unbound = []
    
    for j in range(k_folds):
        if j != i:      # we want to train the model on k-1 other buckets
            training_data_bound.extend(bound_folds[j])      # cmobining all bins from k-1 other buckets to create a "master list" for taining
            training_data_unbound.extend(unbound_folds[j])
            
    b_counts = {}
    u_counts = {}
    
    for nmer in m_plus_1_mers:
        b_counts[nmer] = 1     # pseudocount = 1
        u_counts[nmer] = 1
        
    fill_counts(training_data_bound, b_counts, n)   # getting the frequencies of the m+1 mers using the function defined in step 4
    fill_counts(training_data_unbound, u_counts, n)
    
    b_matrix = transition_prob_matrix(b_counts, m)  # building the transition probability matrices from the training data using function defined in step 5
    u_matrix = transition_prob_matrix(u_counts, m)
    
    # Test on the removed bucket
    test_sequences = bound_folds[i] + unbound_folds[i]      # combinig the 2 removed bins
    test_labels = [1] * len(bound_folds[i]) + [0] * len(unbound_folds[i])   # labeling the bins as '1' for bound and '0' for unbound as we need to know which category it came from for cross-validation
    
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
            
        fold_scores.append(log_p_b - log_p_u)   # calculating the log likelihood score using the transition probability matrices
    
    all_fold_test_labels.append(test_labels)
    all_fold_test_scores.append(fold_scores)

print("Cross-validation scoring is complete")

# Step 7 : Constructing the ROC and Precision-Recall Curves for all 'k' folds

from sklearn.metrics import roc_curve, auc, precision_recall_curve
import matplotlib.pyplot as plt

roc_aucs = []
pr_aucs = []

plt.figure(figsize=(15, 6),dpi=300)

# Subplot 1: ROC Curves
plt.subplot(1, 2, 1)
for i in range(k_folds):
    fpr, tpr, _ = roc_curve(all_fold_test_labels[i], all_fold_test_scores[i])
    fold_auc = auc(fpr, tpr)
    roc_aucs.append(fold_auc)
    plt.plot(fpr, tpr, label=f'Fold {i+1} (AUC = {round(fold_auc, 3)})')

plt.plot([0, 1], [0, 1], 'k--', label='Random Guess')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curves (k-folds)')
plt.legend(loc='lower right') 

# Subplot 2: Precision-Recall Curves
plt.subplot(1, 2, 2)
for i in range(k_folds):
    precision, recall, _ = precision_recall_curve(all_fold_test_labels[i], all_fold_test_scores[i])
    fold_pr_auc = auc(recall, precision)
    pr_aucs.append(fold_pr_auc)
    plt.plot(recall, precision, label=f'Fold {i+1} (AUC = {round(fold_pr_auc, 3)})')

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curves (k-folds)')
plt.legend(loc='upper right') 

plt.tight_layout()
plt.savefig('results.png', bbox_inches='tight', dpi=300)
print("✅ Plots saved as 'results.png' - download to view!")
plt.close()
