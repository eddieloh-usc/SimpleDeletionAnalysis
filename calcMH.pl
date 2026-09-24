#!/usr/bin/perl
# script to calculate microhomology and slippage calls

use strict;

#declare and initialize variables
my $inputfile = $ARGV[0] or die "require input filename\n";
my $fafile = $ARGV[1] or die "require genome reference fasta filename\n";
my $outputfile = $ARGV[2] or die "require output filename\n";

#open the output file
unless (open(OUTFILE, ">$outputfile")){
   print STDERR "Cannot open file \"$outputfile\" to write to!!\n\n";
   exit;
}

### Read the fasta file
print STDERR "reading fasta file...";
my %chrhash=();
my $countchr=0;
my $currchr="";
open (FAFILE,"$fafile") or die "Cannot open reference fasta file!!! \n";
while (my $line=<FAFILE>) {
   chomp $line;
   my $line2=$line;
   if ($line =~ /^>(\S+)/gm) {
      $currchr=$1;
      $countchr++;
      print STDERR ".";
   }else {
      if (exists $chrhash{$currchr}) {
         $chrhash{$currchr}=$chrhash{$currchr}.$line;
      }else {
         $chrhash{$currchr}=$line;
      }
   }
   #last if ($countchr >2);
}
print STDERR "\n";

#printout stored seqname and lengths
foreach my $key (sort {$a cmp $b} (keys %chrhash)) {
   my $seq=$chrhash{$key};
   my $seqlen=length($seq);
   print STDERR "stored\t$key\t$seqlen\n";
}
print STDERR "\n$countchr chrs stored...\n\n";

### Read and process the input file
print STDERR "processing input file...";
my $countline=0;
my $countskipheader=0;
my $countskipcrypt=0;
my $countskipchr=0;
my $countkeep=0;
my %samplehash=();
my %crypthash=();
my %countchrhash=();
open (INPUT,"$inputfile") or die "Cannot open input file!!! \n";
while (my $line=<INPUT>) {
   chomp $line;
   $line =~ s/\r//g;
   my $line2=$line;
   my $line3=$line;
   my $line4=$line;
   my $line5=$line;

   #customize as necessary to input file format
   my @array=split(/\t/,$line);
   my $arraylen=scalar(@array);
   die "numcol error... $arraylen ne 25...\n" if ($arraylen != 25);
   $countline++;
   my $tempchr=$array[2];
   my $tempsample=$array[17];
   my $tempcrypt=$array[18];
   my $tempstart=$array[3]-1; # convert to 0-based coordinates
   my $tempend=$array[9];
   my $tempdellen=$tempend-$tempstart;
   my $tempdelbases=uc($array[12]);
   
   #print header for output file
   if ($tempchr eq "Chromosome") {
      $countskipheader++;
      print OUTFILE "$line2\tdeletelength\trightbases\tmcall\tscall\n";
      #print STDERR "$line2\tdeletelength\trightbases\tmcall\tscall\n";
      next;
   }
   
   #skip non-autosomal chromosomes
   if ((!($tempchr =~ /^chr/gm))||($tempchr eq "chrX")||($tempchr eq "chrY")||($tempchr eq "chrM")) {
      $countskipchr++;
      next;
   }

   $countchrhash{$tempchr}=$tempchr;
   $samplehash{$tempsample}=$tempsample;
   
   #skip unused crypts
   if (($tempcrypt eq "0005_T6")||($tempcrypt eq "0006_T6")||($tempcrypt eq "0012_T2")) {
      $countskipcrypt++;
      next;
   }
   $crypthash{$tempcrypt}=$tempcrypt;
   
   #sanity check on deletion length
   my $tempdellen2=length($tempdelbases);
   die "del length mismatch...$tempdellen...$tempdellen2...\n" if ($tempdellen ne $tempdellen2);

   #extract deletion seq from reference and verify identical to reported deletion seq 
   die "chr $tempchr not found!!!\n" if (!exists($chrhash{$tempchr}));
   my $tempchrseq=$chrhash{$tempchr};
   my $delbases=uc(substr($tempchrseq,$tempstart,$tempend-$tempstart));
   die "extracted variant mismatch...$tempdelbases...$delbases...$countline\n" if ($tempdelbases ne $delbases);
  
   #extract adjacent 50 bases
   my $rightbases=uc(substr($tempchrseq,$tempend,50));

   #assign Mcall
   my @delbases=split("",$delbases);
   my @rightbases=split("",$rightbases);
   my $mcall="";
   for (my $i=1;$i<=$tempdellen;$i++) {
      my $del1=shift(@delbases);
      my $right1=shift(@rightbases);
      if ($del1 eq $right1) {
         $mcall=$mcall."M";
      }else {
         $mcall=$mcall."X";
	 last;
      }
   }

   #assign Scall
   my @delbases=split("",$delbases);
   my @rightbases=split("",$rightbases);
   my $scall="";
   my %uniqdelbase=();
   for (my $i=1;$i<=$tempdellen;$i++) {
      my $del1=shift(@delbases);
      $uniqdelbase{$del1}="$del1";
   }
   if (scalar(keys %uniqdelbase)==1) {
      my $sdel=(keys %uniqdelbase)[0];
      $scall="";
      for (my $i=1;$i<=$tempdellen;$i++) {
         my $del1=$sdel;
         my $right1=shift(@rightbases);
         if ($del1 eq $right1) {
            $scall=$scall."S";
         }else {
            $scall=$scall."X";
	    last;
         }
      }
   }else {
      $scall="NA";
   }

   #print to output 
   print OUTFILE "$line\t$tempdellen\t$rightbases\t$mcall\t$scall\n";
   #print STDERR "$tempchr $tempstart $tempend $tempdelbases $tempdellen $delbases $rightbases $mcall $scall\n";
   $countkeep++;

   print STDERR "\t$countline" if ($countline%1000==0);
   #last if ($countline>=10);
}

#final reporting
my $numsample=scalar(keys %samplehash);
my $numcrypt=scalar(keys %crypthash);
my $numchr=scalar(keys %countchrhash);
print STDERR "\n\n$countline lines read... $countskipheader headerskipped $countskipchr chrskipped... $countskipcrypt cryptskipped... $countkeep processed...\n";
print STDERR "$countkeep processed from $numsample samples, $numchr chrs and $numcrypt crypts...\n"; 

print STDERR "\nComplete!!!\n\n";







