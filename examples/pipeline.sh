#!/bin/bash
RootDir=$(cd `dirname $(readlink -f $0)`; pwd)

### Config start ###
PRIMER3CONFIG="/usr/users/celldev/luf/local/primer3/v2.4.0/x86_64/src/primer3_config/"
### Config ends ###





### Need samtools, muscle, primer3_core

if [ ! -d $PRIMER3CONFIG ]; then
    echo "Error: please modify line5 to specify PRIMER3CONFIG" >&2
    exit 100
fi
homoeologsfile="$RootDir/homoeologs"
Primer3Refseq="$RootDir/ref.fa"
Primer3gff="$RootDir/ref.gff3"

gff3file=/usr/users/celldev/luf/test/scaffolding/4.genesynteny/10.qRT-PCR/all.gff3
NumLines=1;
RunDir=$PWD

samtools faidx $Primer3Refseq
### homoeologsfile: each column are cDNA IDs from a plant/subgenome
###                 each line are homoeologs
### I need the second column seqID to design primers
while read -r nameline; do
	declare -a seqids=()
	seqids=($(echo "$nameline" | perl -lane 'foreach $x (@F) {print $x;}'))
	export NumLines
	echo "Line$NumLines: $nameline"
	cd $RunDir/
	if [ -d $RunDir/${seqids[0]} ]; then
		rm -rf $RunDir/${seqids[0]} >/dev/null 2>&1
	fi
	mkdir -p $RunDir/${seqids[0]}
	cd $RunDir/${seqids[0]}
### Collect all the seqIDs each line to a file
### I renamed the seqIDs as the original IDs are too long to see them in MSF using GeneDoc'
	samtools faidx $Primer3Refseq "${seqids[0]}" > cluster.$NumLines.fa
	samtools faidx $Primer3Refseq "${seqids[1]}" | perl -lne 'if (/^>/) {$newname=">BAC3DL_".$ENV{"NumLines"}; $_=~s/^>/$newname /;} print;' >> cluster.$NumLines.fa
	samtools faidx $Primer3Refseq "${seqids[4]}" | perl -lne 'if (/^>/) {$newname=">EI3DL_".$ENV{"NumLines"}; $_=~s/^>/$newname /;} print;' >> cluster.$NumLines.fa
	samtools faidx $Primer3Refseq "${seqids[3]}" | perl -lne 'if (/^>/) {$newname=">EI3B_".$ENV{"NumLines"}; $_=~s/^>/$newname /;} print;' >> cluster.$NumLines.fa
	samtools faidx $Primer3Refseq "${seqids[2]}" | perl -lne 'if (/^>/) {$newname=">EI3AL_".$ENV{"NumLines"}; $_=~s/^>/$newname /;} print;' >> cluster.$NumLines.fa
### Collect exon length from a GFF3 file, which I used to extract the cDNA sequences
### Be careful this might not work in your case ###
	ExonArr=$(grep "${seqids[1]}" $Primer3gff | perl -lane 'next unless ($F[2] =~/^exon$/i); $strand{$F[6]}++;$len{$F[3]}=$F[4]-$F[3]+1; END {@strandarr=keys %strand; unless (scalar(@strandarr)==1){print "NaN"; exit 0;} if ($strandarr[0] eq "+"){@arr=sort {$a<=>$b} keys %len;}elsif ($strandarr[0] eq "-"){@arr=sort {$b<=>$a} keys %len;}else {print "NaN"; exit 0;} foreach $x (@arr) {push (@arr2, $len{$x});} print join(",", @arr2);}')
	echo "Exonarr: $ExonArr"
### Run HomoeologPrimer.pl
### Col-1, Col2, and Col5 are DD subgenome sequences. So their alleles need to be consistent
### Col-3 and Col4 are AA and BB aleles
### I need to find the DD-specific SNPs, and design primers based on Col-2 sequences
	perl $RootDir/../HomoeologPrimer.pl -i cluster.$NumLines.fa -k ${seqids[0]},BAC3DL_$NumLines,EI3DL_$NumLines -d EI3B_$NumLines,EI3AL_$NumLines -p cluster.$NumLines -m BAC3DL_$NumLines --specific --general --primer3config $PRIMER3CONFIG --exon_len $ExonArr --force > HomoeologPrimer.running.log 2> HomoeologPrimer.running.err
	((NumLines++));
done < $homoeologsfile
