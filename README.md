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
  Yilei Fu - coding  

## Goals
* To develop a tool to detect and genotype (in terms of length) STR in long reads (_de novo_) without the need of genome annotation beforehand
* This tool should be applicable in mammals and plants (basically any phased eukaryotic assembly)

## Introduction/Description
Short tandem repeats (STRs) are motifs of multiple nucleotides in length that are repeated. Due to their repetitive nature, the are subject to high mutation rates and cause several genetic diseases including more than 40 neurological and developmental disorders. Using short-read sequencing data to identify and characterize STRs have been met with some shortcomings including biases introduced by use of PCR. Long-read sequencing can be used to identify their length more accurately than short reads as reads can span across the entire repeat region. Long-reads still have a high error rate.

Although tools have been developed to address the high error-rate problem, they still have limitations such as not being able to consider multiple STRs in a single reads. To adress these issues, we present _STRdust_, a tool to (_de novo_) detect and genotype STRs in long-read sequencing data without prior genome annotation that can be applied in mammals and plants.

## How does it work?  
![alt text](https://raw.githubusercontent.com/collaborativebioinformatics/STRdust/main/STRdustFlowchart.png)

## How to use?  
To run:  
`python STRdust/STRdust.py path_to_bamfile`  

```
usage: Genotype STRs from long reads [-h] [-d DISTANCE] bam

positional arguments:
  bam                   phased bam file

optional arguments:
  -h, --help            show this help message and exit
  -d DISTANCE, --distance DISTANCE
                        distance across which two events should be merged
 ```

## Quickstart

### Input  
  * Phased bam alignment file  

### Output  
  * vcf file with STR genotype calls  

| #chrom | start | end | repeat_seq | size |
| --------- | :-------: | :-----: | :------------: | -------: |
| 22 | 10521893 |	10521909 |	GCA |	16 |
| 22 |  10522628 |	10522647 |	AATA |	19 |
| 22 |	10534566 |	10534586 |	GTTTT |	20 |
| 22 |	10511934 |	10511957 |	AG |	23 |
| 22 |	10516550 |	10516595 |	TA |	45 |
| 22 |	10516550 |	10516572 |	TATATATG |	22 |


## Testing  
This tool was tested using simulated reads for human chromosome 22 and tomato chromosome 1. Long reads were simulated from the GRCh38 (human) and SL4.0 (tomato) reference genome assemblies by first simulating two haplotypes for STRs and furtherthermore with SNPs. The simulated reads were then used to create a phased bam that was used as input for STRdust.   

## Installation  
Installation using a conda environment

1. Clone git repository  
`git clone https://github.com/collaborativebioinformatics/STRdust.git`  

2. Change to git repository folder  
`cd <repository-folder-path>`

3. Install environment  
`conda create --name STRdust -c bioconda python=3.7 pysam spoa pyspoa mreps`  

4. Activate environment  
`conda activate STRdust`  

5. Dry run  
`python STRdust/STRdust.py -h`  


 
