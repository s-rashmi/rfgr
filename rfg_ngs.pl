#!/usr/bin/perl

use Getopt::Long;
use Cwd 'abs_path';

my $input = "";
my $r1 = "";
my $r2 = "";
my $output = "";
my $kmer = 45;
my $pe = 0;

GetOptions(
  "reads=s" => \$input,
  "k=i" => \$kmer,
  "pe!" => \$pe,
  "r1=s" => \$r1,
  "r2=s" => \$r2,
);

if( $r1 ne "" && $r2 ne "" ){
  $pe = 1;
}

if( !$pe ){
  print "Read set:\n$input\n";
} else {
  print "Read set:\n$r1\n";
  print "Read set:\n$r2\n";
}

my $kmerfile = "";
if( $pe ){
  $kmerfile = $r1."_k".$kmer;
} else {
  $kmerfile = $input."_k".$kmer;
}
$uniqkmerfile = $kmerfile."_uniq";
$kmerfreqfile = $kmerfile."_freq";
if( $pe ){
  $nrkmerfile_r1 = $r1."_k".$kmer."_nr";
  $nrkmerfile_r2 = $r2."_k".$kmer."_nr";
} else {
  $nrkmerfile = $kmerfile."_nr";
}

# generate k-mers
print "K-mers of genome in:\n$kmerfile\n"; 
print "Generating k-mers of length $kmer from input....\n";

if( !$pe ){
  open FASTQ, "$input";
} else {
  open FASTQ, "$r1";
}

open OUT, ">$kmerfile";
$line = 0;
$size = 0;
$count = 1;
%seen;
while( <FASTQ> ){
  chomp;
  $line++;
  if( $line == 2 ){
    $size = length($_) - $kmer; 
    for($i = 0; $i < $size; $i++){
      $j = $i+1;
      $sub = substr($_, $i, $kmer);
      $rc = revC($sub);
      print OUT ">".$count."_".$j."_".$kmer."-mer\n";
      if( exists $seen{$rc} ){
        print OUT "$rc\n";
        $seen{$rc}++;
      } else {
        print OUT "$sub\n";  
        $seen{$sub}++;
      }
    } 
  }
  if( $line == 4 ){
    $line = 0;
    $count++;
  }
}
close FASTQ;
if( $pe ){
  open FASTQ, "$r2";
  while( <FASTQ> ){
    chomp;
    $line++;
    if( $line == 2 ){
      $size = length($_) - $kmer; 
      for($i = 0; $i < $size; $i++){
        $j = $i+1;
        $sub = substr($_, $i, $kmer);
        $rc = revC($sub);
        print OUT ">".$count."_".$j."_".$kmer."-mer\n";
        if( exists $seen{$rc} ){
          print OUT "$rc\n";
          $seen{$rc}++;
        } else {
          print OUT "$sub\n";  
          $seen{$sub}++;
        }
      } 
    }
    if( $line == 4 ){
      $line = 0;
      $count++;
    }
  }
  close FASTQ;
} 
close OUT;
print "DONE\n";

# extracting unique k-mers
print "Extracting unique k-mers\n";

open UNIQ, ">$uniqkmerfile";
open FREQ, ">$kmerfreqfile";
foreach $k (keys %seen){
  if( $seen{$k} == 1 ){
    print UNIQ "$k\n";
    delete $seen{$k};
  } else {
    print FREQ "$k\t".$seen{$k}."\n";
  }
}
close FREQ;
close UNIQ;
print "DONE\n";

# get median frequency of non-unique k-mers
$m = median(values %seen);
print "Median Frequency: $m\n";

# tag reads with repeatitive k-mers
# repetitive = 3.25 * median frequency
# Many k-mers are unique because of sequencing error
# thus, to adjust for the median frequency
# we take 3.25 times
%rep;
$thres = 3.25 * $m;
foreach $k (keys %seen){
  if( $seen{$k} >= $thres){
    $rep{$k}++;
  }
}
undef %seen;

if( !$pe ){
  open FASTQ, "$input";
  open OUT, ">$nrkmerfile";
} else {
  open FASTQ, "$r1";
  #open OUT, ">$nrkmerfile_r1";
}
$line = 0;
$size = 0;
$flag = 1;
@print = ();
$ref = {};
$mc = 0;
# total k-mers in a read is L - K + 1
# if more than 1/4 of the k-mers are repetitive
# then the read is tagged as repetitive

while( <FASTQ> ){
  chomp;
  $line++;
  push @print, $_;
  if( $line == 2 ){
    $mc = 0;
    $size = length($_) - $kmer; 
    FOR: for($i = 0; $i < $size; $i++){
      $sub = substr($_, $i, $kmer);
      $rc = revC($sub);
      if( exists $rep{$sub} || exists $rep{$rc} ){
      #  $mc++;
      #}
      #if( $mc >= 0.125 * $size ){
        $flag = 0;
        last FOR;
      }
    } 
  }
  if( $line == 4 ){
    $line = 0;
    if( $flag == 1){
      if( !$pe ){
        print OUT join("\n",@print),"\n";
      } else {
        push @{$ref->{$print[0]}}, @print; 
      }
    }
    $flag = 1;
    @print = ();
  }
}
close FASTQ;
if( !$pe ){
  close OUT;
}

$mc = 0;
$flag = 1;
my $count = 0;
if( $pe ){
  open FASTQ, "$r2";
  open OUT2, ">$nrkmerfile_r2";
  open OUT1, ">$nrkmerfile_r1";
  while( <FASTQ> ){
    chomp;
    $line++;
    push @print, $_;
    if( $line == 2 ){
      $mc = 0;
      $size = length($_) - $kmer; 
      FOR: for($i = 0; $i < $size; $i++){
        $sub = substr($_, $i, $kmer);
        $rc = revC($sub);
        if( exists $rep{$sub} || exists $rep{$rc} ){
        #  $mc++;
        #}
        #if( $mc >= 0.125 * $size ){
          $flag = 0;
          last FOR;
        }
      } 
    }
    if( $line == 4 ){
      $line = 0;
      if( $flag == 1){
        if( exists $ref->{$print[0]} ){
          print OUT1 join("\n",@{$ref->{$print[0]}}),"\n";
          print OUT2 join("\n",@print),"\n";
          $ref->{$print[0]} = ""; 
        } else {
          $count++;
        }
      } else {
        $count++;
      }
      $flag = 1;
      @print = ();
    }
  }
  close FASTQ;
  close OUT1;
  close OUT2;
}

print "$count pairs of repetitive reads were removed\n";

##### Subroutines ######
sub revC {
  my $seq = $_[0];
  my $comp = $seq;
  $comp =~ tr/ATGC/TACG/;
  my $rev = reverse($comp);

  return $rev;
}

sub median {
  my @a = sort {$a <=> $b} @_;
  my $mid = int @a/2;
  return (@a % 2 ) ? $a[$mid] : int ($a[$mid] + $a[$mid -1])/2.0;
}


