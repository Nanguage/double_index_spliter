# Split double indexed fastq file

This program split the fastq file to serveral parts according to it's 
double index sequence.

Double index reads example:
```
@E00477:236:HGGN3CCXY:1:1101:4229:1467 1:N:0:NCAGGCGA+NTCCTTAC
NGCCACGCCATGCCAACCTACACCGGCAGAGCGNCAGCAGTCTAGCTGTGNCAGTTTTGGGGTCCCTGGGGGCCTGAGGGCCGGTAGNCTTACCTGGATACATATTGTTGCTGTCTCTTATACCCCCCTCCGAGCCCACGAGACTAAGGC
+
#A<A-JFFJF7JJJ<FJJJJJ7JJFF--777-7#F777-7-FA-FJF<-7#--7----7AF7AAFJ-7A<JFJJJ<-7FFA-A<7-7#FFFA<JJJJ7A-A-A---77FJJ7F-F-A7--A-A-F--7)AFJ))77J)7A--<F--777A
```

In this reads the index sequence is "NCAGGCGA+NTCCTTAC", the index-a is 'NCAGGCGA' and index-b is 'NTCCTTAC'.

## Dependancy

Python >= 2.7

``` bash
$ pip install click cutadapt biopython
```

## Usage:

```
Usage: double_index_spliter.py [OPTIONS] INPUT_FASTQ INDEX_FILE

Options:
  -m, --mismatch TEXT     mismatch threshold of index-a and index-b, i.e. -m
                          2,3 means 2 bases mismatch on index-a, and 3 bases
                          mismatch on b
  -O, --outdir TEXT       path to output splited fastq files.
  -z, --gzip / --no-gzip  compress output fastq files with gzip.
  --phred INTEGER         encode of fastq quality string.
  --help                  Show this message and exit.

```

example:
``` bash
$ python double_index_spliter.py example/exsample.fq example/example_index.txt -O output/
```

## Index config file

Index config file format example:

```
index-1  GCTACGCT    CTCCTTAC
index-2  CGAGGCTG    CTCCTTAC
index-3  CCAGGCGA    ATCCTTAC
```

It has three columns: index-name, index-a, index-b