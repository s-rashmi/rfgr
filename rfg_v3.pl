#!/usr/bin/perl

use Getopt::Long;
use Cwd 'abs_path';

my $input = "";
my $contig = "";
my $output = "";
my $kmer = 16;
my $gff = "";
my $mm = 0; 

GetOptions(
  "genome=s" => \$input,
  "contig=s" => \$contig,
  "k=i" => \$kmer,
  "gff=s" => \$gff,
  "mm=i" => \$mm
);

my $input_path = "";
if( $input ne "" ){
  $input_path = abs_path($input);
  print "Genome:\n$input_path\n";
} elsif ( $contig ne "" ){
  $input_path = abs_path($contig);
  print "Contigs:\n$input_path\n";
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

$kmerfile = $input_path."_k".$kmer;
$kmerfile1 = $input1."_k".$kmer;
$mapfile = $kmerfile."_map_m".$mm;
$mapfile1 = $kmerfile1."_map_m".$mm;
$logfile = $kmerfile."_log";
$output = $mapfile."_repeats";
$output1 = $mapfile1."_repeats";
$histfile = $mapfile."_hist";

#print "$input\n$input_path\n$input1\n$output\n$output1\n";

$bowtie_path = "/home/rashmi/Documents/Project/bowtie-1.1.2/bowtie-build"; 
@char = split(//, $bowtie_path);
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
$bowtie_path1 = join('',@char);

# build bowtie index
print "Builing bw-index...";
system ("$bowtie_path1 -q $input1 $input1"); 
print "DONE\n";

# generate k-mers
print "Generating $kmer-mers in:\n$kmerfile\n"; 

open FASTA, "$input_path";
$seq = {};
$id = "";
$s = "";
while( $line = <FASTA> ){
  chomp $line;
  if( $line =~ /^>/ ){
    if( $id ne ""){
      $seq->{$id} = $s;
      $id = $line;
      $s = "";
    } else {
      $id = $line;
    }
  } else {
    $line = uc($line);
    $s .= $line;
 }
  # to register the last entry
  $seq->{$id} = $s;
}
close FASTA;

# print length($seq),"\n"; 
# for UTI89 5065741

my $size = "";
my $id = "";
my $name = "";

open OUT, ">$kmerfile";
foreach $id (keys %$seq ){
  # The parentheses provide list context, in which 
  # a pattern match returns the captured parts.  
  # In scalar context, the return value is just 
  # true/false depending on whether the regex matched
  $size = length($seq->{$id}) - $kmer;
  next if $size < $kmer;
  
  if( $input ne "" ){
    my %seen = (); 
    ($name) = $id =~ /(>\w*\|?\d+)/;
    for( $i = 0; $i <= $size; $i++ ){
      $j = $i+1;
      $sub = substr($seq->{$id}, $i, $kmer);
      #$rc = revC($sub);
      #if( !exists $seen{$sub} && !exists $seen{$rc}){
        print OUT $name."_".$j."\n"; 
        print OUT "$sub\n";
      #}
      #$seen{$sub} = 1;
      #$seen{$rc} = 1;
      #print "$i\n" if $i % 10000 == 0;
    }
  }
  if( $contig ne "" ){
    ($name) = $id =~ /(>NODE_\d+)/;
    for( $i = 0; $i <= $size; $i++ ){
      $j = $i+1;
      $sub = substr($seq->{$id}, $i, $kmer);
      print OUT $name."_".$j."\n"; 
      print OUT "$sub\n";
    }
  }
}
close OUT;
print "DONE\n";

print "Mapping repetitive $kmer-mers to reference with atmost $mm mismatch...";
print "\n$mapfile\n"; 
$bowtie_path = "/home/rashmi/Documents/Project/bowtie-1.1.2/bowtie"; 

@char = split(//, $bowtie_path);
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
$bowtie_path1 = join('',@char);

# bowtie run
# bowtie -f -v 0 -a -m 2 -B 1 NC_007646.1.txt NC_007646.1.txt_lmers repeatHits
# Bowtie options
# f = fasta
# v = end-to-end hits with <= v mismatches
# a = all alignments per read
# m = suppress alignments if <m hits (hacked option)
# B = offbase set to 1

$bowtie_output = `$bowtie_path1 -f -v $mm -a -m 2 -B 1 $input1 $kmerfile1 $mapfile1 2>&1`;  
print "DONE\n";
print "Bowtie Log in: \n$logfile\n";

open LOG, ">$logfile";
print LOG $bowtie_output;
close LOG;

my $BTREE;
if( $gff ne ""){ 
  my $gff_path = abs_path($gff);
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
  my $gff1 = join('',@char);

  open GFF, "$gff_path";
  while( $line = <GFF> ){ 
    next if ($line =~ /^#/);
    next if ($line !~ /gene/);
  
    @f = split("\t", $line);
    next if ($f[2] ne "gene");
  
    $g->{$f[3]}->{end} = $f[4];
    $g->{$f[3]}->{strand} = $f[6];
 
    @desc = split(";", $f[8]);
    @n = split("=", $desc[1]);

    if( $n[0] eq "Name" ){
      $g->{$f[3]}->{name} = $n[1];
    } else {
      @n = split("=", $desc[0]);
      $g->{$f[3]}->{name} = $n[1];
    }
  }
  close GFF;
  @genes = keys %$g;
  $BTREE = make_btree($g, 0, scalar(@genes));
}

print "Filtering mapped $kmer-mers...\n";
open MAP, "$mapfile";

$ref = {};
# EDIT: make it a tie hash - to get ordered hash - no need to sort again

$freq = {};
%seen;
$pos1 = "";
$id = "";
$o = "";
$r = "";
$p = "";
while( $line = <MAP>){
  ($id, $strand, $r, $pos) = (split('\t', $line))[0..3]; 
  #id strand ref pos seq quality extra

  if( $input ne "" ){
    @spl = split('_', $id); #gi|556_205  #pos of origin
    $id = $spl[0];
    $ori = $spl[1]; #gi|556_205  #pos of origin
    ($name) = $r =~ /(\w*\|?\d+)/;
  } elsif( $contig ne "") {
    @spl = split('_', $id); #NODE_1_1
    $ori = $spl[2];
    $id = join("_", $spl[0], $spl[1]);
    ($name) = $r =~ /(NODE_\d+)/;
  }

  $flag = 1;
  $pos = -$pos if( $strand eq "-" );
  next if($id eq $name && $ori == $pos);
  if( exists $ref->{$id} ){
    if( exists $ref->{$id}->{$name} ){
      F: for($i=0; $i < scalar(@{$ref->{$id}->{$name}->{pos}}); $i++){
        if( $ori == ${$ref->{$id}->{$name}->{run_ori}}[$i] + 1
            && $pos == ${$ref->{$id}->{$name}->{run_pos}}[$i] +1 ){
          ${$ref->{$id}->{$name}->{k}}[$i] = ${$ref->{$id}->{$name}->{k}}[$i] + 1;
          ${$ref->{$id}->{$name}->{run_ori}}[$i] = ${$ref->{$id}->{$name}->{run_ori}}[$i] + 1;
          ${$ref->{$id}->{$name}->{run_pos}}[$i] = ${$ref->{$id}->{$name}->{run_pos}}[$i] + 1;
          $flag = 0;
          last F;
        }
      }
      if( $flag ){
        push @{$ref->{$id}->{$name}->{ori}}, $ori; 
        push @{$ref->{$id}->{$name}->{run_ori}}, $ori; 
        push @{$ref->{$id}->{$name}->{pos}}, $pos; 
        push @{$ref->{$id}->{$name}->{run_pos}}, $pos; 
        push @{$ref->{$id}->{$name}->{k}}, $kmer; 
      }
    } else {
      push @{$ref->{$id}->{$name}->{ori}}, $ori; 
      push @{$ref->{$id}->{$name}->{run_ori}}, $ori; 
      push @{$ref->{$id}->{$name}->{pos}}, $pos; 
      push @{$ref->{$id}->{$name}->{run_pos}}, $pos; 
      push @{$ref->{$id}->{$name}->{k}}, $kmer; 
    }
  } else {
    push @{$ref->{$id}->{$name}->{ori}}, $ori; 
    push @{$ref->{$id}->{$name}->{run_ori}}, $ori; 
    push @{$ref->{$id}->{$name}->{pos}}, $pos; 
    push @{$ref->{$id}->{$name}->{run_pos}}, $pos; 
    push @{$ref->{$id}->{$name}->{k}}, $kmer; 
  }
}
close MAP;
print "DONE\n";

$count = {};
$c= 0; $d = 0;
open OUT, ">$output";
foreach $id (sort {$a cmp $b} keys %$ref){
  my %seen = ();
  foreach $r (sort {$a cmp $b} keys %{$ref->{$id}}){
    @ori = sort {$a <=> $b} @{$ref->{$id}->{$r}->{ori}};
    for($i=0; $i < scalar(@ori); $i++){
      $type = ""; 
      $ng_o = "";
      $p = ${$ref->{$id}->{$r}->{pos}}[$i];
      $rp = ${$ref->{$id}->{$r}->{run_pos}}[$i];
      $o = ${$ref->{$id}->{$r}->{ori}}[$i];
      $k = ${$ref->{$id}->{$r}->{k}}[$i];
      if( $contig ne "" ){
        if( $id eq $r ){
          next if( exists $seen{$o} );
          $seen{$o} = 1;
          $seen{abs($p)} = 1;
          $c++;
        } else {
          $d++;
        } 
      }
      if( $input ne ""){
        if( $id eq $r ){
          # find closest gene
          if( $gff ne ""){
            if( abs($p) - $o <= 500 || abs($p) + $k - $o <= 500 ){
              $ng_o = closest_gene($BTREE, $o);
            }
          }
  
          # classify repeat
          if( $p < 0 ){
            $type = "Inverted";
            if( abs($rp) == $o ){
              $type .= "- palindrome repeats";
              $count->{inv}->{pal}++;
            } elsif( abs($rp) < $o + $k ){
              $type .= " overlapping repeats";
              $count->{inv}->{ovr}++;
            } else {
              $type .= " dispersed repeats";
              $count->{inv}->{dis}++;
            } 
          } else {
            $type = "Direct";
            if( $p == $o + $k ){
              $type .= "- tandem repeats";
              $count->{dir}->{tan}++;
            } elsif( $p < $o + $k ){
              $type .= " overlapping repeats";
              $count->{dir}->{ovr}++;
            } else {
              $type .= " dispersed repeats";
              $count->{dir}->{dis}++;
            } 
          }
        }
      }
      if( $p < 0 ){
        if( $input ne ""){
          print OUT join("\t", $id, $o."-".($o+$k-1), $k, $r, (abs($p)+$kmer-1)."-".abs($rp), $type, $g->{$ng_o}->{name}),"\n";  
        } else {
          print OUT join("\t", $id, $o."-".($o+$k-1), $k, $r, (abs($p)+$kmer-1)."-".abs($rp)),"\n";  
        }
      } else {
        if( $input ne "" ){
          print OUT join("\t", $id, $o."-".($o+$k-1), $k, $r, $p."-".($p + $k-1), $type, $g->{$ng_o}->{name}),"\n";  
        } else {
          print OUT join("\t", $id, $o."-".($o+$k-1), $k, $r, $p."-".($p + $k-1)),"\n";  
        }
      }
    }
  }
}
close OUT;

if( $input ne ""){
  open HIST, ">$histfile";
  foreach $t (keys %$count){
    $dir = "Direct" if( $t eq "dir");
    $dir = "Inverted" if( $t eq "inv");
  
    foreach $sub (keys %{$count->{$t}}){
      $ty = "Overlapping" if($sub eq "ovr");
      $ty = "Tandem" if($sub eq "tan");
      $ty = "Palindrome" if($sub eq "pal");
      $ty = "Dispersed" if($sub eq "dis");
  
      print HIST join("\t", $dir, $ty, $count->{$t}->{$sub}), "\n";
    }
  }
  close HIST;
} elsif ( $contig ne ""){
  print $c+($d/2)," pairs of repeat locations\n";
} 

#################
## Subroutines ##
#################
sub revC {
  my $seq = $_[0];
  my $comp = $seq;
  $comp =~ tr/ATGC/TACG/;
  my $rev = reverse($comp);

  return $rev;
}

sub make_btree {
  my $h = $_[0];
  my $first = $_[1];
  my $nof = $_[2];
  my @k = sort{ $a <=> $b } keys %$g;

  if( $nof == 1){
    return { Value => $k[$first] },
  } else {
    my $m = int($nof/2);
    my $val = ($k[$first + $m -1] + $k [$first +$m])/2;
    return {
      Value => $val,
      Left => make_btree($h,$first, $m),
      Right => make_btree($h,$first+$m, $nof-$m),
    }
  }
}

sub closest_gene {
  my $btree = $_[0];
  my $key = $_[1];

  while( exists $btree->{Left} ){
    my $v = $btree->{Value};
    $btree = $btree->{ $key < $v ? 'Left' : 'Right' }
  }
  return $btree->{Value};
}
