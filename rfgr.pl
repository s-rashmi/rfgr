#!/usr/bin/perl

use Getopt::Long;
use Cwd 'abs_path';

my $input = "";
my $contig = "";
my $reads = "";
my $r1 = "";
my $r2 = "";
my $pe = 0;
my $output = "";
my $kmer = 16;
my $gff = "";

GetOptions(
  "genome=s" => \$input,
  "contig=s" => \$contig,
  "read=s" => \$reads,
  "r1=s" => \$r1,
  "r2=s" => \$r2,
  "k=i" => \$kmer,
  "gff=s" => \$gff,
);


my $input_path = "";
my $input1 = "";
if( $input ne "" ){
  $input_path = abs_path($input);
  print "Genome:\n$input_path\n";
} elsif ( $contig ne "" ){
  $input_path = abs_path($contig);
  print "Contigs:\n$input_path\n";
} elsif( $reads ne "" ){
  $input_path = abs_path($reads);
}

@char = split(//, $input_path);
for($i=0; $i < scalar(@char); $i++){
  if( $char[$i] =~ /\(/ ){
    $char[$i] = "\"".$char[$i];
  } elsif( $char[$i] =~ /\)/){
    $char[$i] = $char[$i]."\"";
  } elsif( $char[$i] =~ /\n/ ){
    $char[$i] = "";
  } elsif( $char[$i] !~ /[a-zA-Z0-9-_\.]/ && $char[$i] !~ /\//){
    $char[$i] = "\\".$char[$i];
  }
}
$input1 = join('',@char);

my $input_path1 = "";
my $input_path2 = "";
my $input2 = "";
if( $r1 ne "" && $r2 ne "" ){
  $pe = 1;
  $input_path1 = abs_path($r1);

  @char1= split(//, $input_path1);
  for($i=0; $i < scalar(@char1); $i++){
    if( $char1[$i] =~ /\(/ ){
      $char1[$i] = "\"".$char1[$i];
    } elsif( $char1[$i] =~ /\)/){
      $char1[$i] = $char1[$i]."\"";
    } elsif( $char1[$i] =~ /\n/ ){
      $char1[$i] = "";
    } elsif( $char1[$i] !~ /[a-zA-Z0-9-_\.]/ && $char1[$i] !~ /\//){
      $char1[$i] = "\\".$char1[$i];
    }
  }
  $input1 = join('',@char1);

  $input_path2 = abs_path($r2);

  @char2 = split(//, $input_path2);
  for($i=0; $i < scalar(@char2); $i++){
    if( $char2[$i] =~ /\(/ ){
      $char2[$i] = "\"".$char2[$i];
    } elsif( $char2[$i] =~ /\)/){
      $char2[$i] = $char2[$i]."\"";
    } elsif( $char2[$i] =~ /\n/ ){
      $char2[$i] = "";
    } elsif( $char2[$i] !~ /[a-zA-Z0-9-_\.]/ && $char2[$i] !~ /\//){
      $char2[$i] = "\\".$char2[$i];
    }
  }
  $input2 = join('',@char2);
}

my $gff_path = "";
my $gff1 = "";
if( $gff ne ""){
  $gff_path = abs_path($gff);
  @char = split(//, $gff_path);
  for($i=0; $i < scalar(@char); $i++){
    if( $char[$i] =~ /\(/ ){
      $char[$i] = "\"".$char[$i];
    } elsif( $char[$i] =~ /\)/){
      $char[$i] = $char[$i]."\"";
    } elsif( $char[$i] =~ /\n/ ){
      $char[$i] = "";
    } elsif( $char[$i] !~ /[a-zA-Z0-9-_\.]/ && $char[$i] !~ /\//){
      $char[$i] = "\\".$char[$i];
    }
  }
  $gff1 = join('',@char);
}

my $run = "";
my $retcode = "";
if( $input ne "" & $gff ne "" ){
  $run = `perl rfg_suffix.pl -genome $input1 -k $kmer -gff $gff1`;
  $retcode = $?;
  if( $retcode == 0){
    print $run,"\nComplete!\n";  
  } else {
    print $run,"\nError: $retcode\n";  
  }
} elsif( $input ne "" ){
  $run = `perl rfg_suffix.pl -genome $input1 -k $kmer`;
  $retcode = $?;
  if( $retcode == 0){
    print $run,"\nComplete!\n";  
  } else {
    print $run,"\nError: $retcode\n";  
  }
} elsif( $contig ne "" ){
  $run = `perl rfg_v3.pl -contig $input1 -k $kmer`;
  $retcode = $?;
  if( $retcode == 0){
    print $run,"\nComplete!\n";  
  } else {
    print $run,"\nError: $retcode\n";  
  }
} 
if( $pe ){
  if( $kmer == 16 ){
    $kmer = 51;
  }
  $run = `perl rfg_ngs.pl -r1 $input1 -r2 $input2 -k $kmer`;
  $retcode = $?;
  if( $retcode == 0){
    print $run,"\nComplete!\n";  
  } else {
    print $run,"\nError: $retcode\n";  
  }
} elsif( $reads ne "" ){
  if( $kmer == 16 ){
    $kmer = 51;
  }
  $run = `perl rfg_ngs.pl -reads $input_path -k $kmer`;
  $retcode = $?;
  if( $retcode == 0){
    print $run,"\nComplete!\n";  
  } else {
    print $run,"\nError: $retcode\n";  
  }
} elsif( $r1 ne "" && $r2 eq "" ){
  print "ERROR: Only Forward reads provided! Please provide Reverse reads also\n";
  exit;
} elsif( $r1 eq "" && $r2 ne "" ){
  print "ERROR: Only Reverse reads provided! Please provide Forward reads also\n";
  exit;
}
