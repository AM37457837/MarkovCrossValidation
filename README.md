# Markov Model TF Binding Classifier
 This project implements a Markov model–based classifier to predict whether a given transcription factor (TF) binds to a given 200 bp DNA sequence bin, using k-fold cross-validation and evaluation via ROC and precision–recall curves.

## 1. Requirements
 • Language: Python 3.8+
 • Python libraries: pandas, biopython, matplotlib, scikit-learn, gdown
 • Install via: pip install pandas biopython matplotlib scikit-learn gdown

## 2. File Structure
``` text
.
├── MarkovCrossValidation.py     # Main script (run this; includes timestamps)
├── simpler_version.py           # Simplified version
├── fasta_download.py            # Download chrX.fasta from NCBI
├── tsv_download.py              # Download chrX_200bp_bins.tsv from Google Drive
├── output_for_chromosome_4/     # Sample chr4 outputs (ROC_3.png, PR_3.png, chr4_results_summary.txt)
├── results/                     # (Generated: plots)
└── README.md
```
## 3. Data Download
• FASTA (chr1.fasta to chr22.fasta): Run python fasta_download.py, enter chromosome number (1-22).( Downloads from NCBI using RefSeq accessions.)

• TSV (chr1_200bp_bins.tsv to chr22_200bp_bins.tsv): Run python tsv_download.py, enter chromosome number (1-22).( Downloads via gdown from shared Google Drive files.)

• Place files in working directory (or current folder) before running main script. Paths are relative—no absolutes used.

## 4. Input Format
### 4.1 TSV file
Named chr(num)_200bp_bins.tsv. Columns: start , end, CTCF (B/U), REST (B/U), EP300 (B/U).
Example header:
```text
start	end	CTCF	REST	EP300
0	200	B	U	U
```
### 4.2 FASTA file
• Named chr(num).fasta in same directory. 
• Script extracts 200 bp bins using TSV coordinates.

## 5. How to Run
• Run from repo root/terminal (after downloading data):
```text
python MarkovCrossValidation.py
Follow interactive prompts:
TSV filename (e.g., chr1_200bp_bins.tsv).
TF: CTCF, REST, or EP300.
Folds k: 3–5.
Markov order m: 0–10.
```

Example workflow (for chr1, CTCF, k=5, m=6; see chr4 sample in output_for_chromosome_4/):

```text
python tsv_download.py          # Enter 1 → chr1_200bp_bins.tsv
python fasta_download.py        # Enter 1 → chr1.fasta
python MarkovCrossValidation.py # Enter: chr1_200bp_bins.tsv, CTCF, 5, 6
```
## 6. What the Script Does
• Obtains the FASTA sequences of the bins in the tsv file
• Segregates the bins into 'bound' and 'unbound' groups based on the TF of interest
• Randomly divides the bins into 'k' buckets such that each bucket gets roughly equal number of bins
• List which contains all possible 'm+1 mers'
• Counting all the occurences of m+1 mers in the bins
• Building the Transition Probability Matrix
• Training the model on k-1 buckets and using the remaining 1 for cross-validation
• Constructing the ROC and Precision-Recall Curves for all 'k' folds

## 7. Outputs
 --Console: Bin counts, fold progress, total execution time .
 --Plots: ROC_Curves_for_{k}_fold_validation_for_order_{m}.png + Precision-Recall_Curves_for_{k}_fold_validation_for_order_{m} (in png format).
 --Sample: output_for_chromosome_4 (chr4 ROC_3.png, PR_3.png, chr4_results_summary.txt).

## 8. Installation
```text
git clone https://github.com/AM37457837/classifier.git
cd classifier
pip install pandas biopython matplotlib scikit-learn gdown
```
## 9. Notes/Limitations
 Supports CTCF/REST/EP300 only .
 Assumes exact 200 bp bins, matching TSV/FASTA coords.
 Chr 1-22 only; simpler_version.py for quick tests.
