#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use FindBin qw($Bin);
use FuhaoPerl5Lib::AlignKit qw/ReadMsf MsfSnpLocation/;
use FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/;
use Data::Dumper qw /Dumper/;
use File::Basename;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa --keep seqID1,seqID2 --diff seqID3,seqID4 \
        --output mu.out --exon_len INT1,INT2,INT3 \
        primer3config [Options]

Version: v20180622

Requirements:
    Programs[in PATH]: 
        muscle: www.drive5.com/muscle/
        primer3_core: github.com/primer3-org/primer3
    Modules: 
        Scalar::Util, Cwd, Getopt::Long, FindBin
        Data::Dumper, File::Basename
        FuhaoPerl5Lib: github.com/lufuhao/FuhaoPerl5Lib

Descriptions:
    Design primers to amplify specfic homoeolog in wheat
    as well as general primers to amplify all the seqs
    Steps (Included):
        1. Align multiple fasta file using 'muscle'
        2. Identify (--primerseq)-specific SNPs
        3. Use these SNPs as the last base at primer 3'end
        4. Collect primers
    Skills:
        1. do privide primer3config, otherwise primer3 fails
           usually /(Primer3Root)/src/primer3_config
             (Primer3Root) is the folder you installed primer3
        2. use --force to generate some primers even if it 
           violates specific constraints if you can not get 
           any without this option. 
           AND you need to evaluate the primer using some
           programs, like multiple-primer-analyzer from
           thermofisher
        3. Need to specify --exon_len if you want to know
           whether primers cross exons

Options:
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
    --version|v!
        Print current SCRIPT version;

Example:
    perl $0 

Author:
    Fu-Hao Lu
    Post-Doctoral Scientist in Micheal Bevan laboratory
    Cell and Developmental Department, John Innes Centre
    Norwich NR4 7UH, United Kingdom
    E-mail: Fu-Hao.Lu\@jic.ac.uk

EOH
###HELP ends#########################################################
die USAGE unless @ARGV;



###Receving parameter################################################
my ($help, $verbose, $debug, $version);
my ($inputfa);
my $keep_seqs; my $diff_seqs;
my $seqs_config={};
my $outpfx='';
my $output='';
my $primerseq='';
my $test_specific=0;
my $test_general=0;
my $path_to_primer3config;
my $str_exon_len='NaN';
my $seqlen='NaN';
my $specific_target=1;
my $specific_force=0;
my $gap_delimiter='.';### gap_delimiter for alignments, it's dot for muscle
my $max_num_primers=5;
my $amplicon_size='75-250';
my $primer_tm_min=55;
my $primer_tm_max=65;
my $primer_tm_opt=60;
my $primer_len_min=18;
my $primer_len_max=25;
my $primer_len_opt=21;

GetOptions(
	"help|h!" => \$help,
	"input|i=s" => \$inputfa,
	"output|o:s" => \$output,
	"keep|k=s" => \$keep_seqs,
	"diff|d=s" => \$diff_seqs,
	"prefix|p:s" => \$outpfx,
	"primerseq|m=s" => \$primerseq,
	"specific!" => \$test_specific,
	"general!" => \$test_general,
	'primer3config:s' => \$path_to_primer3config,
	'exon_len:s' => \$str_exon_len,
	'force' => \$specific_force,
	'numprimers|n' => \$max_num_primers,
	'size' => \$amplicon_size,
	'tmopt' => \$primer_tm_opt,
	'tmmin' => \$primer_tm_min,
	'tmmax' => \$primer_tm_max,
	'lenopt' => \$primer_len_opt,
	'lenmin' => \$primer_len_min,
	'lenmax' => \$primer_len_max,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);
$verbose=0 unless (defined $verbose);

unless (defined $path_to_primer3config) {
	if (exists $ENV{'PRIMER_THERMODYNAMIC_PARAMETERS_PATH'}) {
		$path_to_primer3config=$ENV{'PRIMER_THERMODYNAMIC_PARAMETERS_PATH'}
	}
}
unless ($path_to_primer3config=~/^\/.*\/$/ and -d $path_to_primer3config and -s $path_to_primer3config."/stack.ds") {
	die "Error: invalid primer3config path: $path_to_primer3config\n";
}
my @exons=split (',', $str_exon_len);
foreach my $indnum (@exons) {
	unless (defined $indnum and $indnum=~/^\d+$/ and $indnum>0) {
		print STDERR "Error: invalid Exon length array: $str_exon_len\n";
		@exons=(); last;
	}
}
my %exonborder=();
my $exonstart=1;
my $exonend=0;
foreach my $indnum (@exons) {
	$exonend+=$indnum;
	$exonborder{$exonstart}{'min'}=$exonstart;
	$exonborder{$exonstart}{'max'}=$exonend;
	$exonstart+=$indnum;
}




### input and output ################################################
$inputfa=abs_path($inputfa);
unless (defined $inputfa and -s $inputfa) {
	die "Error: invalid input fasta file\n";
}
my @arr1=split(/,/, $keep_seqs);
foreach (@arr1) {
	${$seqs_config}{'keep'}{$_}++;
}
my @arr2=split(/,/, $diff_seqs);
foreach (@arr2) {
	${$seqs_config}{'diff'}{$_}++;
}
if ($outpfx eq '') {
	$outpfx=basename($inputfa);
}
if ($output eq '') {
	$output=$outpfx.".primers"
}
unless (defined $primerseq and $primerseq =~/^\S+$/) {
	die "Error: please specify primer template --primerseq\n";
}
if ($specific_force) {
	$specific_target=0;
}
unless (defined $max_num_primers and $max_num_primers=~/^\d+$/ and $max_num_primers >0) {
	die "Error: invalid number of primers: --numprimers|-n\n";
}
unless (defined $amplicon_size and $amplicon_size=~/^\d+-\d+$/) {
	die "Error: invalid Amplicon size[70-250]\n";
}
unless (defined $primer_tm_opt and $primer_tm_opt=~/^\d+$/) {
	die "Error: invalid optimum Tm value: --tmopt\n";
}
unless (defined $primer_tm_min and $primer_tm_min=~/^\d+$/) {
	die "Error: invalid minimum Tm value: --tmmin\n";
}
unless (defined $primer_tm_max and $primer_tm_max=~/^\d+$/) {
	die "Error: invalid maximum Tm value: --tmmax\n";
}
unless (defined $primer_len_opt and $primer_len_opt=~/^\d+$/) {
	die "Error: invalid optimum primer length value: --lenmin\n";
}
unless (defined $primer_len_min and $primer_len_min=~/^\d+$/) {
	die "Error: invalid minimum primer length value: --lenmin\n";
}
unless (defined $primer_len_max and $primer_len_max=~/^\d+$/) {
	die "Error: invalid maximum primer length value: --lenmax\n";
}

if ($debug) {
	print "### Main: Debug \$seqs_config\n";
	print Dumper $seqs_config;
	print "\n";
}


### Main ############################################################
print "##### SUMMARY #####\n";
print "Primer3_config : $path_to_primer3config\n";
print "Input   Fasta  : $inputfa\n";
print "Output  prefix : $outpfx\n";
print "Output primers : $output\n\n";
print "Seq to keep    : ", join(",", @arr1), "\n";
print "Seq to diff    : ", join(",", @arr2), "\n\n";
print "Primer template: ", $primerseq, "\n";
print "Exon length    : ", join (",", @exons), "\n";
print "Primer specific:   : ". ($test_specific==1) ? "Yes" : "NO" . "\n";
if ($test_specific) {
	print "Number primers (specific)  : ", $max_num_primers, "\n";
}
print "Primer general :   : ". ($test_general==1) ? "Yes" : "NO" . "\n";
if ($test_general) {
	print "Number primers (general)  : ", $max_num_primers, "\n";
}
print "\n";
print "Amplicon size  : ", $amplicon_size, "\n";
print "Primer Tm      : MIN $primer_tm_min OPT $primer_tm_opt MAX $primer_tm_max\n";
print "Primer length  : MIN $primer_len_min OPT $primer_len_opt MAX $primer_len_max\n";

if ($debug) {
	print "### Main: Debug \%exonborder\n";
	print Dumper \%exonborder;
	print "\n";
}
print "\n";
print "\n";

print "Main: Info: 1. run muscle\n" if ($verbose);
my $msfout="$outpfx.msf";
unless (exec_cmd_return("muscle -msf -in $inputfa -out $msfout -quiet")) {
	die "Error01: muscle running failed\n";
}

print "Main: Info: 2.Read MSF\n" if ($verbose);
my ($test1, $seqalgn)=ReadMsf($msfout);
unless ($test1) {
	die "Error02: ReadMsf running failed\n";
}
unless (exists ${$seqalgn}{$primerseq}) {
	die "Error: Primer template sequence should be one of the aligned sequences\n";
}
if ($debug) {
	print "### Main: Debug \$seqalgn\n";
	print Dumper $seqalgn;
	print "\n";
}

print "Main: Info: 3. Get specific and shared alleles\n" if ($verbose);
my ($test2, $pos2allele, $sharedallele)=MsfSnpLocation($seqalgn, $seqs_config);
unless ($test2) {
	die "Error03: MsfSnpLocation running failed\n";
}
if ($debug) {
	print "### Main: Debug : \$pos2allele\n";
	print Dumper $pos2allele;
	print "\n";
}
if ($debug) {
	print "### Main: Debug : \$sharedallele\n";
	print Dumper $sharedallele;
	print "\n";
}


unless (open FINALPRIMER, ">", $output) {
	die "Error: can not write --output\n";
}

if ($test_specific) {
	$seqlen='NaN';
	print FINALPRIMER "################################\n";
	print FINALPRIMER "########### SPECIFIC ###########\n";
	print FINALPRIMER "################################\n";
	print "### Info: $primerseq\n";
	print STDERR "### Info: $primerseq\n";
	unless (scalar(keys %{$pos2allele})>0) {
		die "Error04: no specfic allele detected\n";
	}
	my ($test3, $physicloci, $template)=&MsfLoci2PhysicalLocation(${$seqalgn}{$primerseq}, $pos2allele, $gap_delimiter);
	unless ($test3) {
		die "Error: MsfLoci2PhysicalLocation failed\n";
	}
	foreach my $x (sort {$a<=>$b} keys %{$physicloci}) {
		print FINALPRIMER "### POS $x\n";
		if ($specific_target==1) {
			print FINALPRIMER "### POS $x MODE pick_discriminative_primers\n";
			my $primer3input="$outpfx.specific.pos$x";
			close SPECIFIC if (defined fileno(SPECIFIC));
			unless (open SPECIFIC, ">", $primer3input) {
				die "Error: can not write primer3 [specific] config file: $primer3input\n";
			}
			print SPECIFIC "SEQUENCE_ID=$primerseq\n";
			print SPECIFIC "SEQUENCE_TEMPLATE=$template\n";
			print SPECIFIC "PRIMER_TASK=pick_discriminative_primers\n";
			print SPECIFIC "PRIMER_PICK_LEFT_PRIMER=1\n";
			print SPECIFIC "PRIMER_PICK_INTERNAL_OLIGO=0\n";
			print SPECIFIC "PRIMER_PICK_RIGHT_PRIMER=1\n";
			print SPECIFIC "PRIMER_OPT_SIZE=21\n";
			print SPECIFIC "PRIMER_MIN_SIZE=18\n";
			print SPECIFIC "PRIMER_MAX_SIZE=25\n";
			print SPECIFIC "PRIMER_OPT_TM=$primer_tm_opt\n";
			print SPECIFIC "PRIMER_MAX_TM=$primer_tm_max\n";
			print SPECIFIC "PRIMER_MIN_TM=$primer_tm_min\n";
			print SPECIFIC "PRIMER_MAX_NS_ACCEPTED=0\n";
			print SPECIFIC "PRIMER_PRODUCT_SIZE_RANGE=$amplicon_size\n";
			print SPECIFIC "P3_FILE_FLAG=1\n";
			print SPECIFIC "PRIMER_EXPLAIN_FLAG=1\n";
			print SPECIFIC "PRIMER_NUM_RETURN=$max_num_primers\n";
			print SPECIFIC "SEQUENCE_TARGET=", ($x-1), ",1\n";
			print SPECIFIC "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=$path_to_primer3config\n";
			print SPECIFIC "PRIMER_PICK_ANYWAY=1\n";
			print SPECIFIC "=\n";
			close SPECIFIC;
			unless (-s $primer3input) {
				die "Error: failed ti write primer3 config file\n";
			}
			unless (exec_cmd_return("primer3_core --format_output --output $primer3input.out --error $primer3input.err < $primer3input")) {
				die "Error01: primer3_core running failed\n";
			}
			unless (-s "$primer3input.err") {
				unlink "$primer3input.err" if (-e "$primer3input.err");
			}
			my $test5=&ProcessPrimer3Out("$primer3input.out");
			unlink "$primer3input.out" if ($test5==2);
			unless ($test5) {
				print STDERR "Warnings: get primer3 out error: $primer3input.out\n"
			}
		}
		if ($specific_force==1) {
			print FINALPRIMER "### POS $x MODE force LEFT\n";
			my $primer3input="$outpfx.specific.pos$x.left";
			close SPECIFIC if (defined fileno(SPECIFIC));
			unless (open SPECIFIC, ">", $primer3input) {
				die "Error: can not write primer3 [specific] config file: $primer3input\n";
			}
			print SPECIFIC "SEQUENCE_ID=$primerseq\n";
			print SPECIFIC "SEQUENCE_TEMPLATE=$template\n";
			print SPECIFIC "PRIMER_TASK=generic\n";
			print SPECIFIC "PRIMER_PICK_LEFT_PRIMER=1\n";
			print SPECIFIC "PRIMER_PICK_INTERNAL_OLIGO=0\n";
			print SPECIFIC "PRIMER_PICK_RIGHT_PRIMER=1\n";
			print SPECIFIC "PRIMER_OPT_SIZE=21\n";
			print SPECIFIC "PRIMER_MIN_SIZE=18\n";
			print SPECIFIC "PRIMER_MAX_SIZE=25\n";
			print SPECIFIC "PRIMER_OPT_TM=$primer_tm_opt\n";
			print SPECIFIC "PRIMER_MAX_TM=$primer_tm_max\n";
			print SPECIFIC "PRIMER_MIN_TM=$primer_tm_min\n";
			print SPECIFIC "PRIMER_MAX_NS_ACCEPTED=0\n";
			print SPECIFIC "PRIMER_PRODUCT_SIZE_RANGE=$amplicon_size\n";
			print SPECIFIC "P3_FILE_FLAG=1\n";
			print SPECIFIC "PRIMER_EXPLAIN_FLAG=1\n";
			print SPECIFIC "PRIMER_NUM_RETURN=$max_num_primers\n";
			print SPECIFIC "SEQUENCE_FORCE_LEFT_END=", ($x-1), "\n";
			print SPECIFIC "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=$path_to_primer3config\n";
			print SPECIFIC "PRIMER_PICK_ANYWAY=1\n";
			print SPECIFIC "=\n";
			close SPECIFIC;
			unless (-s $primer3input) {
				die "Error: failed to write primer3 config file\n";
			}
			unless (exec_cmd_return("primer3_core --format_output --output $primer3input.out --error $primer3input.err < $primer3input")) {
				die "Error01: primer3_core running failed\n";
			}
			unless (-s "$primer3input.err") {
				unlink "$primer3input.err" if (-e "$primer3input.err");
			}
			my $test5=&ProcessPrimer3Out("$primer3input.out");
			unlink "$primer3input.out" if ($test5==2);
			unless ($test5) {
				print STDERR "Warnings: get primer3 out error: $primer3input.out\n"
			}
			
			print FINALPRIMER "### POS $x MODE force RIGHT\n";
			$primer3input="$outpfx.specific.pos$x.right";
			close SPECIFIC if (defined fileno(SPECIFIC));
			unless (open SPECIFIC, ">", $primer3input) {
				die "Error: can not write primer3 [specific] config file: $primer3input\n";
			}
			print SPECIFIC "SEQUENCE_ID=$primerseq\n";
			print SPECIFIC "SEQUENCE_TEMPLATE=$template\n";
			print SPECIFIC "PRIMER_TASK=generic\n";
			print SPECIFIC "PRIMER_PICK_LEFT_PRIMER=1\n";
			print SPECIFIC "PRIMER_PICK_INTERNAL_OLIGO=0\n";
			print SPECIFIC "PRIMER_PICK_RIGHT_PRIMER=1\n";
			print SPECIFIC "PRIMER_OPT_SIZE=21\n";
			print SPECIFIC "PRIMER_MIN_SIZE=18\n";
			print SPECIFIC "PRIMER_MAX_SIZE=25\n";
			print SPECIFIC "PRIMER_OPT_TM=$primer_tm_opt\n";
			print SPECIFIC "PRIMER_MAX_TM=$primer_tm_max\n";
			print SPECIFIC "PRIMER_MIN_TM=$primer_tm_min\n";
			print SPECIFIC "PRIMER_MAX_NS_ACCEPTED=0\n";
			print SPECIFIC "PRIMER_PRODUCT_SIZE_RANGE=$amplicon_size\n";
			print SPECIFIC "P3_FILE_FLAG=1\n";
			print SPECIFIC "PRIMER_EXPLAIN_FLAG=1\n";
			print SPECIFIC "PRIMER_NUM_RETURN=$max_num_primers\n";
			print SPECIFIC "SEQUENCE_FORCE_RIGHT_END=", ($x-1), "\n";
			print SPECIFIC "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=$path_to_primer3config\n";
			print SPECIFIC "PRIMER_PICK_ANYWAY=1\n";
			print SPECIFIC "=\n";
			close SPECIFIC;
			unless (-s $primer3input) {
				die "Error: failed to write primer3 config file\n";
			}
			unless (exec_cmd_return("primer3_core --format_output --output $primer3input.out --error $primer3input.err < $primer3input")) {
				die "Error01: primer3_core running failed\n";
			}
			unless (-s "$primer3input.err") {
				unlink "$primer3input.err" if (-e "$primer3input.err");
			}
			my $test6=&ProcessPrimer3Out("$primer3input.out");
			unlink "$primer3input.out" if ($test5==2);
			unless ($test6) {
				print STDERR "Warnings: get primer3 out error: $primer3input.out\n"
			}
		}
	}
}




if ($test_general) {
	$seqlen='NaN';
	print FINALPRIMER "################################\n";
	print FINALPRIMER "########### GENERAL ###########\n";
	print FINALPRIMER "################################\n";
	my ($test4, $physicExclude, $template)=&MsfLoci2ExcludedLocation(${$seqalgn}{$primerseq}, $sharedallele, $gap_delimiter);
	unless ($test4) {
		die "Error: MsfLoci2ExcludedLocation failed\n";
	}
	if ($debug) {
		print "### Main: Debug : \$physicExclude\n";
		print Dumper $physicExclude;
		print "\n";
	}
	my @sequence_target=();
	foreach my $x (keys %{$physicExclude}) {
		if ($x=~/^(\d+)-(\d+)$/) {
			push (@sequence_target, ($1-1).','.($2-$1+1));
		}
		else {
			die "Error: excluded region format error: $x\n";
		}
	}
	my $primer3input="$outpfx.general";
	close SPECIFIC if (defined fileno(SPECIFIC));
	unless (open SPECIFIC, ">", $primer3input) {
		die "Error: can not write primer3 [specific] config file: $primer3input\n";
	}
	print SPECIFIC "SEQUENCE_ID=$primerseq\n";
	print SPECIFIC "SEQUENCE_TEMPLATE=$template\n";
	print SPECIFIC "PRIMER_TASK=generic\n";
	print SPECIFIC "PRIMER_PICK_LEFT_PRIMER=1\n";
	print SPECIFIC "PRIMER_PICK_INTERNAL_OLIGO=0\n";
	print SPECIFIC "PRIMER_PICK_RIGHT_PRIMER=1\n";
	print SPECIFIC "PRIMER_OPT_SIZE=21\n";
	print SPECIFIC "PRIMER_MIN_SIZE=18\n";
	print SPECIFIC "PRIMER_MAX_SIZE=25\n";
	print SPECIFIC "PRIMER_OPT_TM=$primer_tm_opt\n";
	print SPECIFIC "PRIMER_MAX_TM=$primer_tm_max\n";
	print SPECIFIC "PRIMER_MIN_TM=$primer_tm_min\n";
	print SPECIFIC "PRIMER_MAX_NS_ACCEPTED=0\n";
	print SPECIFIC "PRIMER_PRODUCT_SIZE_RANGE=$amplicon_size\n";
	print SPECIFIC "P3_FILE_FLAG=1\n";
	print SPECIFIC "PRIMER_EXPLAIN_FLAG=1\n";
	print SPECIFIC "PRIMER_NUM_RETURN=", 8 * $max_num_primers, "\n";
	print SPECIFIC "SEQUENCE_EXCLUDED_REGION=", join(' ', @sequence_target), "\n";
	print SPECIFIC "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=$path_to_primer3config\n";
	print SPECIFIC "=\n";
	close SPECIFIC;
	unless (-s $primer3input) {
		die "Error: failed to write primer3 config file\n";
	}
	unless (exec_cmd_return("primer3_core --format_output --output $primer3input.out --error $primer3input.err < $primer3input")) {
		die "Error01: primer3_core running failed\n";
	}
	unless (-s "$primer3input.err") {
		unlink "$primer3input.err" if (-e "$primer3input.err");
	}
	my $test5=&ProcessPrimer3Out("$primer3input.out");
	unlink "$primer3input.out" if ($test5==2);
	unless ($test5) {
		print STDERR "Warnings: get primer3 out error: $primer3input.out\n"
	}
}
close FINALPRIMER;


################### sub modules ############################
### Global: $pos2allele, $primerseq,$gap_delimiter, $seqlen
sub MsfLoci2PhysicalLocation {
	my ($MLPLaligned_seq, $MLPLspecific_alleles, $MLPLdelimiter)=@_;
	
	my $MLPLsubinfo='SUB(MsfLoci2PhysicalLocation)';
	my $MLPLloc=0;
	my %MLPLreturn=();
	my $MLPLret_seq='';
	$seqlen=length($MLPLaligned_seq);
	
	unless ($seqlen>0) {
		print STDERR $MLPLsubinfo, "Error: seq length 0\n";
		return 0;
	}
	
	MLPLLOOP: for (my $MLPLx=0; $MLPLx<$seqlen; $MLPLx++) {
		my $MLPLthis_allele=substr($MLPLaligned_seq, $MLPLx, 1);
		if ($MLPLthis_allele eq $MLPLdelimiter) {
			next MLPLLOOP;
		}
		else {
			$MLPLloc++;
			$MLPLret_seq.=$MLPLthis_allele;
			unless ($MLPLthis_allele=~/^[atcgnATCGN]{1,1}$/) {
				print STDERR $MLPLsubinfo, "Warnings: not standard base at position: Align ", $MLPLx+1, " Physical $MLPLloc\n";
			}
		}
		if (exists ${$MLPLspecific_alleles}{$MLPLx+1}) {
			print $MLPLsubinfo, "Info: Align ", $MLPLx+1, " Physical $MLPLloc\n";
			$MLPLreturn{$MLPLloc}=${$MLPLspecific_alleles}{$MLPLx+1}{'keep'};
		}
	}

	return (1, \%MLPLreturn, $MLPLret_seq);
}

sub MsfLoci2ExcludedLocation {
	my ($MLELaligned_seq, $MLELconserved_alleles, $MLELdelimiter)=@_;
	
	my $MLELsubinfo='SUB(MsfLoci2ExcludedLocation)';
	my $MLELloc=0;
	my %MLELreturn=();
	my %MLELexcluded=();
	my $MLELret_seq='';
	$seqlen=length($MLELaligned_seq);
	
	unless ($seqlen>0) {
		print STDERR $MLELsubinfo, "Error: seq length 0\n";
		return 0;
	}
	
	MLELLOOP1: for (my $MLELx=0; $MLELx<$seqlen; $MLELx++) {
		my $MLELthis_allele=substr($MLELaligned_seq, $MLELx, 1);
		if ($MLELthis_allele eq $MLELdelimiter) {
			next MLELLOOP1;
		}
		else {
			$MLELloc++;
			$MLELret_seq.=$MLELthis_allele;
			unless ($MLELthis_allele=~/^[atcgnATCGN]{1,1}$/) {
				print STDERR "Warnings: not standard base at position: Align ", $MLELx+1, " Physical $MLELloc\n";
			}
		}
		
		unless (exists ${$MLELconserved_alleles}{$MLELx+1}) {
			print $MLELsubinfo, "Info: SEQ: $primerseq Align ", $MLELx+1, " Physical $MLELloc\n";
			$MLELexcluded{$MLELloc}++;
		}

	}
	my @MLELarr=sort {$a<=>$b} keys %MLELexcluded;
	my $MLELlast=0;
	my $MLELstart=0;
	MLELLOOP2: for (my $MLELy=0; $MLELy<scalar(@MLELarr); $MLELy++) {
		if ($MLELy==0) {
			$MLELlast=$MLELarr[$MLELy];
			$MLELstart=$MLELarr[$MLELy];
		}
		else {
			unless ($MLELarr[$MLELy] == ($MLELlast+1)) {
				$MLELreturn{"$MLELstart-$MLELlast"}=$MLELlast-$MLELstart+1;
				$MLELstart=$MLELarr[$MLELy];
			}
			$MLELlast=$MLELarr[$MLELy];
		}
	}
	$MLELreturn{"$MLELstart-$MLELlast"}=$MLELlast-$MLELstart+1;
	
	return (1, \%MLELreturn, $MLELret_seq);
}



### Global: %exonborder
sub ProcessPrimer3Out {
	my $PPOprimer3out=shift;
	
	my $PPOsubinfo="SUB(ProcessPrimer3Out)";
	my $PPOlinenum=0;
	my $PPOprimer_suffix=0;
	my %PPOprimers=();
	local *PPOINPUT;
	
	unless (defined $PPOprimer3out and -s $PPOprimer3out) {
		print $PPOsubinfo, "Error: invalid primter3 output\n";
		return 0;
	}
	
	close PPOINPUT if (defined fileno(PPOINPUT));
	unless (open PPOINPUT, "<", $PPOprimer3out) {
		print $PPOsubinfo, "Error: can not open primter3 output\n";
		return 0;
	}
	while (my $PPOline=<PPOINPUT>) {
		$PPOlinenum++;
		chomp $PPOline;
		if ($PPOline=~/^NO\s+PRIMERS\s+FOUND$/) {
			print $PPOsubinfo, "Info: No primers: $PPOprimer3out\n";
			return 2;
		}
#OLIGO            start  len      tm     gc%  any_th  3'_th hairpin seq
#LEFT PRIMER        346   21   60.06   52.38    0.00   0.00    0.00 TCCTGCTCATCATCCTCCAGA
#RIGHT PRIMER       474   21   59.93   52.38    0.00   0.00   44.66 CCTTCACAACAGATGGTGGGA
		if ($PPOline=~/LEFT PRIMER\s+\d+/) {
			$PPOprimer_suffix++;
			$PPOline=~s/^.*LEFT PRIMER\s+//;
			my @PPOarr=split(/\s+/, $PPOline);
			$PPOprimers{$PPOprimer_suffix}{'left'}{'note'}=join("\t", @PPOarr);
			$PPOprimers{$PPOprimer_suffix}{'left'}{'pos'}=$PPOarr[0]+$PPOarr[1]-1;
		}
		elsif ($PPOline=~/RIGHT PRIMER\s+\d+/) {
			$PPOline=~s/^.*RIGHT PRIMER\s+//;
			my @PPOarr=split(/\s+/, $PPOline);
			$PPOprimers{$PPOprimer_suffix}{'right'}{'note'}=join("\t", @PPOarr);
			$PPOprimers{$PPOprimer_suffix}{'right'}{'pos'}=$PPOarr[0]-$PPOarr[1]+1;
		}
	}
	close PPOINPUT;
	my @PPOborder=();
	foreach my $PPOa (sort {$a<=>$b} keys %exonborder) {
		push (@PPOborder, $exonborder{$PPOa}{'min'}.'-'.$exonborder{$PPOa}{'max'});
	}
	print FINALPRIMER "### Exons: ", join(",", @PPOborder), "\n";
	print FINALPRIMER "SeqID\tSeqLen\tIndex\tF/R\tstart\tlen\ttm\tgc%\tany_th\t3'_th\thairpin\tseq\tAcross_Exons\n";
	foreach my $PPOx (sort {$a<=>$b } keys %PPOprimers) {
		my %PPOexonid=();
		my $PPOacross_exons=0;
		if (exists $PPOprimers{$PPOx}{'left'} and exists $PPOprimers{$PPOx}{'right'}) {
			foreach my $PPOy (sort {$a<=>$b} keys %exonborder) {
				if ($PPOprimers{$PPOx}{'left'}{'pos'}>=$exonborder{$PPOy}{'min'} and $PPOprimers{$PPOx}{'left'}{'pos'}<=$exonborder{$PPOy}{'max'}) {
					$PPOexonid{$PPOy}++;
				}
				if ($PPOprimers{$PPOx}{'right'}{'pos'}>=$exonborder{$PPOy}{'min'} and $PPOprimers{$PPOx}{'right'}{'pos'}<=$exonborder{$PPOy}{'max'}) {
					$PPOexonid{$PPOy}++;
				}
			}
			if (scalar(keys %PPOexonid)>1) {
				$PPOprimers{$PPOx}{'left'}{'note'}.="\tACROSS_EXONS";
				$PPOprimers{$PPOx}{'right'}{'note'}.="\tACROSS_EXONS";
			}
			else {
				$PPOprimers{$PPOx}{'left'}{'note'}.="\tNO";
				$PPOprimers{$PPOx}{'right'}{'note'}.="\tNO";
			}
			print FINALPRIMER $primerseq, "\t", $seqlen, "\t", $PPOx, "\tF\t", $PPOprimers{$PPOx}{'left'}{'note'}, "\n";
			print FINALPRIMER $primerseq, "\t", $seqlen, "\t", $PPOx, "\tR\t", $PPOprimers{$PPOx}{'right'}{'note'}, "\n";
		}
		elsif (exists $PPOprimers{$PPOx}{$PPOx}{'left'}) {
			print FINALPRIMER $primerseq, "\t", $seqlen, "\t", $PPOx, "\tF\t", $PPOprimers{$PPOx}{'left'}{'note'}, "\tNaN\n";
		}
		elsif (exists $PPOprimers{$PPOx}{'right'}) {
			print FINALPRIMER $primerseq, "\t", $seqlen, "\t", $PPOx, "\tR\t", $PPOprimers{$PPOx}{'right'}{'note'}, "\tNaN\n";
		}
	}
	
	return 1;
}
