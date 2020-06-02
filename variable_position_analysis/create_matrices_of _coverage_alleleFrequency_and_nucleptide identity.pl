use File::Find;

#This script scans directories of VCF and COVERAGE files and yields matrices of coverage, allele frequency and nucleotide identity
#for the examined samples. Each row in the matrices contains at least one sample which shows deviation of a given nucleotide in specified genomic location
#compared to the reference genome

#This script was written by Hila Gingold

#INPUT: VCF and Coverage Files
#OUTPUT: 3 tables:
#Column #1 = chromosome
#Column #2 = coordinate
#Column #3 and on = samples IDs

#nuc_outfile: each cell denotes the identity of the nucleotide in sample X for a given genomic location (as defined in column #1 & #2); "-" for position X which is not covered in sample Y
#cov_outfile: each cell denotes the coverage of the nucleotide in sample X for a given genomic location (as defined in column #1 & #2)
#af_outfile: each cell denotes the allele frequency of detected SNP in sample X for a given genomic location (as defined in column #1 & #2); "NaN" for position X which shows no deviation from the reference genome in sample Y

#OPTIONAL_DESCRIPTION refers to alignment type and other properties

my $cov_outfile='cov_at_pos_OPTIONAL_DESCRIPTION.csv';
my $af_outfile='af_of_pos_OPTIONAL_DESCRIPTION.csv';
my $nuc_outfile='nuc_of_pos_OPTIONAL_DESCRIPTION.csv';

open OUTFILE1, ">$cov_outfile" or die "Can't open $cov_outfile for write:$!\n";

open OUTFILE2, ">$af_outfile" or die "Can't open $af_outfile for write:$!\n";

open OUTFILE3, ">$nuc_outfile" or die "Can't open $nuc_outfile for write:$!\n";


my %allSNP = ();
my %allCOV = ();
my %snpREF = ();

my @sample_list = ();

my @files_vcf = glob("./vcf/*.vcf");   #directory containing VCF file per sample

foreach my $filename (@files_vcf)
{
   open (INFILE_vcf,$filename) or die "Can't open file type $!\n"; 
   if  ($filename =~ m/.\/vcf\/(.*?)\.vcf/)  { 
      my $sample = $1;
      print STDOUT $sample."\n";
      push @sample_list, $sample;  
      while (<INFILE_vcf>) {
         if ($_ =~ /^(.*?)\t+(.*?)\t+(.*?)\t+(.*?)\t+(.*?)\t+100\t+PASS\t+AF=(.*?);(.*?);COV=(.*?)\n/)  {
            my $chrom = $1;
            my $pos = $2;
            my $cov = $8;
            my $af = $6;
            my $ref = $4;
            my $altr = $5;
            my $id = $chrom.",".$pos;
            my $isTransversion = 1;
            if ((($ref eq "A") && ($altr eq "G")) || (($ref eq "G") && ($altr eq "A")) || (($ref eq "C") && ($altr eq "T")) || (($ref eq "T") && ($altr eq "C"))) {
               $isTransversion = 0
            } 
            if ($isTransversion) {  #Ignore transitions				
               if (not exists $snpREF{$id}) {
                  $snpREF{$id} = $ref;
               }
               $allSNP{$id}{$sample} = $cov.",".$af.",".$altr; 
            }
         }
       }
      close INFILE_vcf;
   }
}

my @files_perbase = glob("./cov_w_alt/*.tsv");   #directory containing COVERAGE file per sample

foreach my $filename_pb (@files_perbase) {
   if  ($filename_pb =~ m/.\/cov_w_alt\/(.*?)\.tsv/)  {   
      my $sample = $1;
      print STDOUT $sample."\n";
      open (INFILE,$filename_pb) or die "Can't open file type $!\n";  
      while (<INFILE>) {
         if ($_ =~ /^(.*?)\t+(.*?)\t+(.*?)\t+(.*?)\t+0\t+(.*?)\n/)  {
            my $chrom = $1;
            my $pos = $2;
            my $cov_num = $5;
            my $id = $chrom.",".$pos;
            if (exists $allSNP{$id}) {
               $allCOV{$id}{$sample} = $cov_num;
            }
         }
       }
      close INFILE;
   }
}


# printing


print OUTFILE1 "chr,pos";
print OUTFILE2 "chr,pos";
print OUTFILE3 "chr,pos";

foreach $_ (@sample_list) {
    print OUTFILE1 ",".$_;
    print OUTFILE2 ",".$_;
    print OUTFILE3 ",".$_;
}
print OUTFILE1 "\n";
print OUTFILE2 "\n";
print OUTFILE3 "\n";

for my $snp ( sort keys %allSNP ) {
    print OUTFILE1 $snp;
    print OUTFILE2 $snp;
    print OUTFILE3 $snp;
    foreach $_ (@sample_list) {
        my $samp = $_;
        if (exists $allSNP{$snp}{$samp}) {
            my @snp_tag = split /,/,$allSNP{$snp}{$samp};
            print OUTFILE1 ",".$snp_tag[0];
            print OUTFILE2 ",".$snp_tag[1];
            print OUTFILE3 ",".$snp_tag[2];
        } elsif (exists $allCOV{$snp}{$samp}) {
            print OUTFILE1 ",".$allCOV{$snp}{$samp};
            print OUTFILE2 ",NaN";
            print OUTFILE3 ",".$snpREF{$snp};
        } else {
            print OUTFILE1 ",0";
            print OUTFILE2 ",NaN";
            print OUTFILE3 ",-";
        }
    }
    print OUTFILE1 "\n";
    print OUTFILE2 "\n";
    print OUTFILE3 "\n";
}

