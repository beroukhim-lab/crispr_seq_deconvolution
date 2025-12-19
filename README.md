# crispr_seq_deconvolution
This repo describes an analysis pipeline for local analysis of CRISPR-Seq data using crispresso.


### Run the entire pipeline on a directory with fastq files
```
bash run_pipeline.sh \
  -F /path/to/dir/with/fastq/ \
  -a /path/to/table/with/primer/sequences/sample_annotation.txt \
  -o /path/to/output/dir \
  -d pinellolab/crispresso2 \
  -p linux/amd64
```
