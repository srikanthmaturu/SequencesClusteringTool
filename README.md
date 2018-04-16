# SequencesClusteringTool
A sequences similarity tool that uses cosine similarity approach. The tool currently supports both database quering and clustering.

## Software requirements
* cmake 3.4
* make
* git
* gcc >5.4
* OpenMP

## Sub repositories as dependencies
* [Edlib](https://github.com/srikanthmaturu/edlib.git)
* [FALCONN](https://github.com/srikanthmaturu/FALCONN.git)
* [seqan](https://github.com/seqan/seqan.git)

## Instructions to build the tool
```
git clone https://github.com/srikanthmaturu/SequencesClusteringTool.git
cd SequencesClusteringTool
echo "Initialize and clone sub repositories"
git submodule init
git submodule update
mkdir build
cd build
cmake ..
echo "Build program for PI threshold 90"
make seq_anlyzer_DT_1_NL_3_LT_2_NHT_16_NHB_14_NP_16_TH_50_PT_SparseVectorFloat
echo "Build program for PI threshold 80"
make seq_anlyzer_DT_1_NL_3_LT_2_NHT_16_NHB_14_NP_16_TH_60_PT_SparseVectorFloat
echo "Build program for PI threshold 70"
make seq_anlyzer_DT_1_NL_3_LT_2_NHT_32_NHB_14_NP_300_TH_70_PT_SparseVectorFloat
```

## Sample data
Complete SWISSPROT dataset available at [link](ftp://ftp.ncbi.nih.gov/blast/db/FASTA/) <br/>
**Database protein sequences:** `data/swissprot_10000_sequences` <br/>
**Query protein sequences:** `data/swissprot_10000_sequences_5000_sequences` <br/>
**Raw results:** `data/raw_results` <br/>
**tabular results:** `data/tabular_results` (excell) <br/>

## Database querying

Usage:
 ```
 ./executable 1 [database_file] [file_type] [query_file] [file_type] [dataset_type] [data_type] [minPercentIdentity] [parallel]

```

### Command options

* **database_file** - Input protein database sequences
* **file_type** - *fasta* for fasta files <br/>
                *kmers* for kmer files
* **query_file** - Input protein query sequences
* **file_type** - *fasta* for fasta files <br/>
                *kmers* for kmer files
* **dataset_type** - *0*
* **data_type** - *0* for DNA sequences <br/>
                *1* for protein sequences
* **minimum percent identity threshold** - *xx* Ex: 90, 80, 70
* **parallel** - '0' for single thread <br/>
               '1' for multiple threads

### Output results

The above command produces raw results file with name *query_file*_falconn_results.txt

The results should be interpreted as follows:

```
>Sequence: [query_sequence_id]-[length_of_query_sequence]
[match_database_sequence1_id]-[match_database_sequence1_length]-[pi_between_query_&_match_database_sequence1],[match_database_sequence2_id]-[match_database_sequence2_length]-[pi_between_query_&_match_database_sequence2]
```

Ex:
```
>Sequence: 0-559
0-559-100,1-563-93
```
In the example, it shows that query sequence with id 0 match with database sequences with ids 0 & 1.

### Sample results

Sample raw results for input database: `data/swissprot_10000_sequences` and queries: `data/swissprot_10000_sequences_5000_sequences` <br/>
are at `data/raw_data/PI_TH_90/swissprot_10000_sequences_5000_sequences_falconn_results.txt`

## Database clustering

Usage:
```
./executable 0 [sequences_file] [file_type] [dataset_type] [data_type] [minPercentIdentity] [typeOfAlignment] [stepSize]
```

### Command options

* **sequences_file** - Input protein database sequences
* **file_type** - *fasta* for fasta files <br/>
                *kmers* for kmer files
* **dataset_type** - *0*
* **data_type** - *0* for DNA sequences <br/>
                *1* for protein sequences
* **minimum percent identity threshold** - *xx* Ex: 90, 80, 70
* **typeOfAlignment** - *0* for global alignment <br/>
                      *1* for infix alignment
* **stepSize** - *10000* determines how often the underlying index used for clustering needs to be refreshed

### Output results

The above command produces raw results file with name *sequences_file*.clusters.txt


The results should be interpreted as follows:
```
>ClusterId: [cluster_index]
[representative_sequence_id] *
[sequence1_id] at [pi_between_repr_sequence_&_sequence1]%
[sequence2_id] at [pi_between_repr_sequence_&_sequence2]%
```

Ex:
```
>ClusterId: 4795
9413 *
9415 at 93%
9414 at 93%
```

In the example, it shows that cluster with index 4795 contains sequences with ids 9413, 9415, 9414 and sequence 9413 is the cluster representative

### Sample results

Sample raw results for input `data/swissprot_10000_sequences` are at `data/raw_data/PI_TH_90/swissprot_10000_sequences.clusters.txt`

