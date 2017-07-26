#!/bin/bash

# Download TARGET data from UCSC Xena as processed by TOIL Recompute
# https://xenabrowser.net/datapages/?host=https://toil.xenahubs.net

# Gene Expression
wget https://toil.xenahubs.net/download/target_RSEM_Hugo_norm_count.gz
gunzip target_RSEM_Hugo_norm_count.gz

# Phenotype Data
wget https://toil.xenahubs.net/download/TARGET_phenotype.gz
gunzip TARGET_phenotype.gz

