# HomoeologPrimer

## SYNOPSIS:

perl $0 --input my.fa --keep seqID1,seqID2 --diff seqID3,seqID4 --output mu.out --exon_len INT1,INT2,INT3 --primer3config /path/to/primer3config/ [Options]

Version: v20180622



## Requirements:

+ Programs[in PATH]: 

 - muscle:  [www.drive5.com/muscle/](http://www.drive5.com/muscle/)

 -  primer3_core:  [github.com/primer3-org/primer3](https://github.com/primer3-org/primer3)

+ Modules: 
 
 - Scalar::Util, Cwd, Getopt::Long, FindBin, Data::Dumper, File::Basename

 - **FuhaoPerl5Lib**: [github.com/lufuhao/FuhaoPerl5Lib](https://github.com/lufuhao/FuhaoPerl5Lib)

## Descriptions:

### Design primers to amplify

- [x]  specific homoeolog
- [x]  all the homoeologs

### Steps (Included):

1. Align multiple fasta file using 'muscle'

2. Identify `--primerseq`-specific SNPs

3. Use these SNPs as the last base at primer 3'end

4. Run primer3

5. Collect primers

### Skills:

1. Do privide `--primer3config`, otherwise primer3 fails; 

  > usually /(Primer3Root)/src/primer3_config

  > (Primer3Root) is the folder you installed primer3

2. use `--force` to generate some primers even if it violates specific constraints if you can not get any without this option. if you can not get any without this option; And you **need to evaluate the primer** using some programs: like [multiple-primer-analyzer](https://www.googleadservices.com/pagead/aclk?sa=L&ai=DChcSEwjz3rLZsOfbAhUI4RsKHSelC4YYABAAGgJ3bA&ohost=www.google.co.uk&cid=CAESEeD2UmkcERPFfsH_BJZOycHO&sig=AOD64_1jHeP4MI3gCB23HRWy5rxYO1bFgA&q=&ved=0ahUKEwi21qzZsOfbAhWLcRQKHZxLAnEQ0QwIJw&adurl=) (thermofisher)

3. Need to specify `--exon_len` if you want to know whether primers cross exons

4. See examples/README for more primer selection skills 

## Options:

```
--help|-h
        Print this help/usage;
--input|-i  my.fa
       cDNA multi-fasta file for each homoeolog group
--keep|-k  [seqID1,seqID2,...]
        SeqIDs having consistent alleles, comma-delimited
--diff|-d  [seqID3,seqID4,...]
        SeqIDs having different alleles, comma-delimited
--primerseq|-m <seqID>
        Sequence used as template to design primers
--specific
        Design specific primers for (--primerseq|-m)
--general
        Design general primers for all aligned sequences
--output|-o  <output>
        Final primer output file
--exon_len  <exon1len,exon2len,...>
        Exon length to analyze if primer across exons
        Comma-delimited INTs; in N->C (cDNA) order
--primer3config <path>
        /path/to/primer3config
--force
        Force to output some primers even if it violates
        specific constraints.
--numprimers|-n <INT[5]>
        Number of SNP specific primers for each SNP
        AND number of general primer would be 8 times
--size <INT-INT[75-250]>
        Amplicon size range
--tmopt <INT[60]>
        Optimum Tm value
--tmmin <INT[55]>
        Minimum Tm value
--tmmax <INT[65]>
        Maximum Tm value
--lenopt <INT[21]>
        Optimum primer length value
--lenmin <INT[18]>
        Minimum primer length value
--lenmax <INT[25]>
        Maximum primer length value
--verbose
        Detailed output for trouble-shooting;
--version|v
        Print current SCRIPT version;
```

## Author:

>Fu-Hao Lu

>Post-Doctoral Scientist in Micheal Bevan laboratory

>Cell and Developmental Department, John Innes Centre

>Norwich NR4 7UH, United Kingdom

>E-mail: Fu-Hao.Lu@jic.ac.uk
