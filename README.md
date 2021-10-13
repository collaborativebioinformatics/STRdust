![alt text](https://raw.githubusercontent.com/collaborativebioinformatics/STRdust/main/STRdust-logo.jpg)

## Contributors

  Anneri Lötter, Guangyi Chen and Susanne P. Pfeifer - Writers  
  Luis Paulin - coding  
  Damaris Lattimer - coding  
  Deepak Choubey  
  Kyuil Cho - Sysadmin  
  Kimberley Billingsley  
  Pavel Avdeyev - coding, dataset generation  
  Simone Cree  
  Wouter De Coster - Team Lead  
  Yilei Fu  

## Goals
* To develop a tool to detect and genotype (in terms of length) STR in long reads (_de novo_) without the need of genome annotation beforehand
* This tool should be applicable in mammals and plants (basically any phased eukaryotic assembly)

## Introduction/Description
Short tandem repeats (STRs) are motifs of multiple nucleotides in length that are repeated. Due to their repetitive nature, the are subject to high mutation rates and cause several genetic diseases including more than 40 neurological and developmental disorders. Using short-read sequencing data to identify and characterize STRs have been met with some shortcomings including biases introduced by use of PCR. Long-read sequencing can be used to identify their length more accurately than short reads as reads can span across the entire repeat region. Long-reads still have a high error rate.

Although tools have been developed to address the high error-rate problem, they still have limitations such as not being able to consider multiple STRs in a single reads. To adress these issues, we present _STRdust_, a tool to (_de novo_) detect and genotype STRs in long-read sequencing data without prior genome annotation that can be applied in mammals and plants.

## How does it work?

![alt text](https://raw.githubusercontent.com/collaborativebioinformatics/STR_Integration/main/Flow%20chart%20group2.jpg)

## How to use?
`conda activate STRdust`  
`python STRdust/STRdust.py <path_to_bamfile>`


## Quickstart

### Input  
  * Phased reference genome  
  * fastq files of long-read sequencing run  


### Output  
  * vcf file with STR genotype calls  

## Testing

## Installation  
Installation with conda environment:  
1. `conda create --name STRdust -c bioconda python=3.7 pysam spoa pyspoa mreps`  
2. `conda activate STRdust`  
3. `python STRdust/STRdust.py`  


 
