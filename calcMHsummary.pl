#!/usr/bin/perl
# script to summarize microhomology and slippage calls

use strict;

#declare and initialize variables
my $inputfile = $ARGV[0] or die "require input filename\n";
my $outputfile = $ARGV[1] or die "require output filename\n";

#open the output file
unless (open(OUTFILE, ">$outputfile")){
   print STDERR "Cannot open file \"$outputfile\" to write to!!\n\n";
   exit;
}
print OUTFILE "delete_length\ttype\tcount\n";

### Read the input file
my $countline=0;
my %dellenhash=();
my %mhash=();
my %shash=();
open (INPUT,"$inputfile") or die "Cannot open input file 1!!! \n";
while (my $line=<INPUT>) {
   chomp $line;
   $line =~ s/\r//g;
   my $line2=$line;
   my $line3=$line;
   my $line4=$line;
   my $line5=$line;

   my @array=split(/\t/,$line);
   my $arraylen=scalar(@array);
   die "numcol error... $arraylen ne 29...\n" if ($arraylen != 29);
   $countline++;
   my $tempchr=$array[2];
   if ($tempchr eq "Chromosome") {
      next;
   }
   my $del=$array[12];
   my $dellen=length($del);
   my $mcall=$array[27];
   my $scall=$array[28];

   if (!exists $dellenhash{$dellen}) {
      $dellenhash{$dellen}=1;
      $mhash{$dellen}=$mcall;
      $shash{$dellen}=$scall;
   }else {
      $dellenhash{$dellen}=$dellenhash{$dellen}+1;
      $mhash{$dellen}=$mhash{$dellen}."::$mcall";
      $shash{$dellen}=$shash{$dellen}."::$scall";
   }

   if ($countline%1000==0) {
      print STDERR "$countline:\t$del $dellen $mcall $scall\n";
   }
   #last if ($countline>=1000);
}
print STDERR "\nComplete!!! $countline lines parsed...\n\n";
sleep(1);

my $countkey=0;
foreach my $key (sort {$a <=> $b} (keys %dellenhash)) {
   $countkey++;
   my @mcalls=split(/::/,$mhash{$key});
   my @scalls=split(/::/,$shash{$key});
   my $dellen=$key;

   #next if ($dellen != 9);

   my %tempm=();
   foreach my $mcall (@mcalls) {
      my $mname=$mcall;
      if (!exists $tempm{$mname}) {
         $tempm{$mname}=1;
      }else {
         $tempm{$mname}=$tempm{$mname}+1;
      }
   }

   my %temps=();
   foreach my $scall (@scalls) {
      my $sname=$scall;
      if (!exists $temps{$sname}) {
         $temps{$sname}=1;
      }else {
         $temps{$sname}=$temps{$sname}+1;
      }
   }

   my %tempcomb=();
   foreach my $key (sort {$a cmp $b}(keys %tempm)) {
      my $combname="";
      if ($key =~ /X/gm) {
	 my $combkey=$key;
	 $combkey=~s/X//gm;
         my $comblen=length($combkey);
         $combname=$comblen."MH";
      }else {
         $combname="Slip";
      }
      $tempcomb{$combname}=$tempm{$key};
   }
   foreach my $skey (sort {$a cmp $b}(keys %temps)) {
      next if ($skey eq "NA");
      next if ($skey eq "X");
      next if (!($skey =~ /X/gm));
      my $svalue=$temps{$skey};
      my $mkey=$skey;
      $mkey=~s/S/M/gm;
      $mkey=~s/X//gm;
      my $mlen=length($mkey);
      my $combname=$mlen."MH";
      if ($svalue>0) {
         $tempcomb{"Slip"}=$tempcomb{"Slip"}+$svalue; 
         $tempcomb{$combname}=$tempcomb{$combname}-$svalue; 
      }
   }

   print STDERR "dellens $key\n$dellenhash{$key}\n@mcalls\n@scalls\n";
   foreach my $mkey (sort {$a cmp $b}(keys %tempm)) {
      print STDERR "   $mkey\t$tempm{$mkey}\n";
   }
   print STDERR "-------------\n";
   foreach my $skey (sort {$a cmp $b}(keys %temps)) {
      print STDERR "   $skey\t$temps{$skey}\n";
   }
   print STDERR "-------------\n";
   foreach my $combkey (sort {$b cmp $a}(keys %tempcomb)) {
      print STDERR "   $combkey\t$tempcomb{$combkey}\n";
      print OUTFILE "$dellen\t$combkey\t$tempcomb{$combkey}\n";
   }
   print STDERR "\n";

   #last if ($countkey==10);
}






