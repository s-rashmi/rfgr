#!/usr/bin/perl

use Tree::Suffix;
use Getopt::Long;

my $input_path = "";
my $output = "";
my $kmer = 16;
my $gff_path = "";

GetOptions(
  "genome=s" => \$input_path,
  "k=i" => \$kmer,
  "gff=s" => \$gff_path,
);

$kmerfile = $input_path."_k".$kmer;
$mapfile = $kmerfile."_map";
$output = $mapfile."_repeats";
$histfile = $mapfile."_hist";
$distfile = $mapfile."_dist";

# reads sequences to hash
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
    $s .= $line;
  }
  # to register the last entry
  $seq->{$id} = $s;
}
close FASTA;

my $BTREE;
if( $gff_path ne ""){ 
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

# initialize suffix tree
# set allow_duplicates to FALSE 
$tree = Tree::Suffix->new;
$bool = 0;
$tree->allow_duplicates($bool);

# To maintain order of input sequences
# for the suffix tree with the output 
# sequence id of the search function
print "Building Suffix Tree....\n";
@ids = keys %$seq;
@seqs = ();
foreach $id ( @ids ){
  push @seqs, $seq->{$id};
}
$tree = Tree::Suffix->new(@seqs);
print "DONE\n";

@prtF;
@prtR;
@tmpF;
@tmpR;
%seenF = ();
%seenR = ();
$geneF = "";
$geneR = "";
$countF = {};
$countR = {};

print "Locating maximal repeats...\n";
open DIST, ">$distfile";
open MAP, ">$mapfile";
for( $t = 0; $t < scalar(@ids); $t++ ){
  # extract name
  ($source) = $ids[$t] =~ /(\w*\|?\d+)/;

  # size of genome
  $size = length($seq->{$ids[$t]}) - $kmer;

  # generate k-mers
  for( $i = 0; $i <= $size; $i++ ){
    # pos of origin in genome
    $ori = $i+1;

    # forward k-mer
    $sub = substr($seq->{$ids[$t]}, $i, $kmer);
    @pos = ($tree->search($sub));
    $s = scalar(@pos);

    # if only self hit for current position, then no increments
    # for existing positions needed, so print them
    if( $s == 1 && scalar(@prtF) != 0 ){
      for($l = 0; $l < scalar(@prtF); $l=$l+5){
        ($o) = (split("-",$prtF[$l+1]))[0];
        ($p1,$p2) = (split("-",$prtF[$l+4]));

        # find closest gene
        print DIST join("\t", abs($p1 - $o), abs($p2-$o)),"\n";
        if( abs($p1 - $o) <= 500 || abs($p2 - $o) <= 500 ){
          if( $gff_path ne ""){
            $ng_pos = closest_gene($BTREE, $o);
            $geneF = $g->{$ng_pos}->{name};
          }
        }

        # classify repeats
        $type = "Direct ";
        if( $p1 == $o + $prtF[$l+2] ){
          $type .= "- tandem repeats"; 
          $countF->{tan}++;
        } elsif( $p1 < $o + $prtF[$l+2] ){
          $type .= "overlapping repeats"; 
          $countF->{ovr}++;
        } else {
          $type .= "dispersed repeats"; 
          $countF->{dis}++;
        }

        print MAP join("\t", @prtF[$l..$l+4], $type, $geneF),"\n";
        $type = "";
        $geneF = "";
      }
      @prtF = ();
      @tmpF = ();
      %seenF = ();
    }

    # if pos has multiple hits, only then need to process
    if( $s > 1 ){
      foreach $p (@pos){ 
        $flag = 1;
        # each element is a ref to an array [str_id, start, end]
        # str_id is the sequence index in the array we input
        # start and end are 0-offset positions

        # ignore if match pos same as pos of origin
        next if $ori-1 == $$p[1];
        # ignore if match pos is less than pos of origin
        next if $ori-1 > $$p[1];

        # increment start and end to fit 1-offset pos
        $$p[1]++;
        $$p[2]++;
        # extract name of match sequence and update @$p
        ($r) = $ids[$$p[0]] =~ /(\w*\|?\d+)/;

        #process positions
        if( scalar(@prtF) == 0 ){
        # if prtF is empty, ie first match
          push @prtF, ($source, $ori."-".($ori+$kmer-1), $kmer, $r, $$p[1]."-".$$p[2]);
          push @tmpF, ($source, $ori, $kmer, $r, $$p[1]);
          $seenF{$$p[1]} = 0;
        } else {
          $pos = $$p[1];

          # if hit is an increment of a previously seen hit
          if( exists $seenF{$pos-1} ){
            $in = $seenF{$pos-1};
             
            if( $ori == $tmpF[$in+1]+1 && $pos == $tmpF[$in+4]+1
                && $r == $tmpF[$in+3] ){
              $flag = 0;

              # update the hash with the increased pos no.
              delete $seenF{$pos-1};
              $seenF{$pos} = $in;

              # increment kmer value
              $prtF[$in+2]++;

              # change end point of seq of origin
              # update the tmp array with unchanged value
              $tmpF[$in+1] = $ori;
              $o = $ori+$kmer-1; 
              $prtF[$in+1] =~ s/\-\d+/\-$o/; 

              # change end point of match seq 
              # update the tmp array with unchanged value
              $tmpF[$in+4] = $pos;
              $prtF[$in+4] =~ s/\-\d+/\-$$p[2]/; 
              #$prtF[$in+4] = ($$p[1]-1)."-".$$p[2]; 
            }
          } else {
            J: for($j=0; $j < scalar(@prtF); $j=$j+5){
              if ( $ori == $tmpF[$j+1] && $r eq $prtF[$j+3]
                   && ($pos != $tmpF[$j+4]+1 || $pos != $tmpF[$j+4])){
              # if hit is for the same k-mer but at different pos

                $flag = 0;
                $seenF{$pos} = scalar(@prtF);
                push @prtF, ($source, $ori."-".($ori+$kmer-1), $kmer, $r, $pos."-".($pos+$kmer-1));
                push @tmpF, ($source, $ori, $kmer, $r, $pos);
                last J;
              } elsif( $ori == $tmpF[$j+1]+1 && $r eq $prtF[$j+3] ){
              # if hit is for the next k-mer but at a non-incrementive pos
                $flag = 0;
                $seenF{$pos} = scalar(@prtF);
                push @prtF, ($source, $ori."-".($ori+$kmer-1), $kmer, $r, $pos."-".($pos+$kmer-1));
                push @tmpF, ($source, $ori, $kmer, $r, $pos);
                last J;
              }
            }
          }
          if( $flag ){
            # if not an increment, print the earlier match
            # empty both array and fill with new match
            for($l = 0; $l < scalar(@prtF); $l=$l+5){
              ($o) = (split("-",$prtF[$l+1]))[0];
              ($p1,$p2) = (split("-",$prtF[$l+4]));
      
              # find closest gene
              print DIST join("\t", abs($p1 - $o), abs($p2-$o)),"\n";
              if( abs($p1 - $o) <= 500 || abs($p2 - $o) <= 500 ){
                if( $gff_path ne ""){
                  $ng_pos = closest_gene($BTREE, $o);
                  $geneF = $g->{$ng_pos}->{name};
                }
              }
      
              # classify repeats
              $type = "Direct ";
              if( $p1 == $o + $prtF[$l+2] ){
                $type .= "- tandem repeats"; 
                $countF->{tan}++;
              } elsif( $p1 < $o + $prtF[$l+2] ){
                $type .= "overlapping repeats"; 
                $countF->{ovr}++;
              } else {
                $type .= "dispersed repeats"; 
                $countF->{dis}++;
              }
              print MAP join("\t", @prtF[$l..$l+4], $type, $geneF),"\n";
              $type = "";
              $geneF = "";
            }
            @prtF = ();
            @tmpF = ();
            %seenF = ();

            $seenF{$pos} = 0;
            push @prtF, ($source, $ori."-".($ori+$kmer-1), $kmer, $r, $pos."-".($pos+$kmer-1));
            push @tmpF, ($source, $ori, $kmer, $r, $pos);
          }
        }
      }
    }

    # reverse complement of kmer
    $rc = revC($sub);
    @pos = ($tree->search($rc));

    # if only self hit for current position, then no increments
    # for existing positions needed, so print them
    if( scalar(@pos) == 0 && scalar(@prtR) != 0 ){
      for($l = 0; $l < scalar(@prtR); $l=$l+5){
        ($o) = (split("-",$prtR[$l+1]))[0];
        ($p1,$p2) = (split("-",$prtR[$l+4]));

        # find closest gene
        print DIST join("\t", abs($p1 - $o), abs($p2-$o)),"\n";
        if( abs($p1 - $o) <= 500 || abs($p2 - $o) <= 500 ){
          if( $gff_path ne ""){
            $ng_pos = closest_gene($BTREE, $o);
            $geneR = $g->{$ng_pos}->{name};
          }
        }

        # classify repeats
        $type = "Inverted ";
        if( $p1 == $o + $prtR[$l+2]-1 ){
          $type .= "- palindrome repeats"; 
          $countR->{pal}++;
        } elsif( $p1 == $o + $prtF[$l+2] ){
          $type .= "- tandem repeats"; 
          $countR->{tan}++;
        } elsif( $p2 < $o + $prtR[$l+2] ){
          $type .= "overlapping repeats"; 
          $countR->{ovr}++;
        } else {
          $type .= "dispersed repeats"; 
          $countR->{dis}++;
        }
        print MAP join("\t", @prtR[$l..$l+4], $type, $geneR),"\n";
        $type = "";
        $geneR = "";
      }
      @prtR = ();
      @tmpR = ();
      %seenR = ();
    } 

    foreach $p (@pos){
      # each element is a ref to an array [str_id, start, end]
      # str_id is the sequence index in the array we input
      # start and end are 0-offset positions

      $flagR = 1;
      # ignore if match pos is less than pos of origin
      # and is diff is less than k-mer --> palindrome
      next if ($ori-1 > $$p[1]) && (($ori-1)-$$p[1] > $kmer);

      # increment start and end to fit 1-offset pos
      $$p[1]++;
      $$p[2]++;
      # extract name of match sequence and update @$p
      ($r) = $ids[$$p[0]] =~ /(\w*\|?\d+)/;

      # process positions
      if( scalar(@prtR)  == 0 ){
      # if prtR is empty, ie first match
        push @prtR, ($source, $ori."-".($ori+$kmer-1), $kmer, $r, $$p[2]."-".$$p[1]);
        push @tmpR, ($source, $ori, $kmer, $r, $$p[1]);
        $seenR{$$p[1]} = 0;
      } else {
        $pos = $$p[1];

        # if hit is an increment of a previously seen hit
        if( exists $seenR{$pos+1} ){
          $in = $seenR{$pos+1};

          if( $ori == $tmpR[$in+1]+1 && $pos == $tmpR[$in+4]-1
              && $r == $tmpR[$in+3] ){
            $flagR = 0;

            # update the hash with the increased pos no.
            delete $seenR{$pos+1};
            $seenR{$pos} = $in;

            # increment kmer value
            $prtR[$in+2]++;

            # change end point of seq of origin
            # update the tmp array with unchanged value
            $tmpR[$in+1] = $ori;
            $o = $ori+$kmer-1;
            $prtR[$in+1] =~ s/\-\d+/\-$o/;

            # change end point of match seq 
            # update the tmp array with unchanged value
            $tmpR[$in+4] = $pos;
            $prtR[$in+4] =~ s/-\d+/\-$pos/;
            #$prtR[$in+4] = $$p[1]."-".($$p[2]+1);
          }
        } else {
          J: for($j=0; $j < scalar(@prtR); $j=$j+5){
            if ( $ori == $tmpR[$j+1] && $r eq $prtR[$j+3]
                 && ($pos != $tmpR[$j+4]-1 || $pos != $tmpR[$j+4])){
            # if hit is for the same k-mer but at different pos

              $flagR = 0;
              $seenR{$pos} = scalar(@prtR);
              push @prtR, ($source, $ori."-".($ori+$kmer-1), $kmer, $r, ($pos+$kmer-1)."-".$pos);
              push @tmpR, ($source, $ori, $kmer, $r, $pos);
              last J;
            } elsif( $ori == $tmpR[$j+1]+1 && $r eq $prtR[$j+3] ){
            # if hit is for the next k-mer but at a non-incrementive pos
              $flagR = 0;
              $seenR{$pos} = scalar(@prtR);
              push @prtR, ($source, $ori."-".($ori+$kmer-1), $kmer, $r, ($pos+$kmer-1)."-".$pos);
              push @tmpR, ($source, $ori, $kmer, $r, $pos);
              last J;
            }
          }
        }
        if( $flagR ){
          # if not an increment, print the earlier match
          # empty both array and fill with new match
          for($l = 0; $l < scalar(@prtR); $l=$l+5){
            ($o) = (split("-",$prtR[$l+1]))[0];
            ($p1,$p2) = (split("-",$prtR[$l+4]));
    
            # find closest gene
            print DIST join("\t", abs($p1 - $o), abs($p2-$o)),"\n";
            if( abs($p1 - $o) <= 500 || abs($p2 - $o) <= 500 ){
              if( $gff_path ne ""){
                $ng_pos = closest_gene($BTREE, $o);
                $geneR = $g->{$ng_pos}->{name};
              }
            }
    
            # classify repeats
            $type = "Inverted ";
            if( $p1 == $o + $prtR[$l+2] ){
              $type .= "- palindrome repeats"; 
              $countR->{pal}++;
            } elsif( $p1 == $o + $prtF[$l+2] ){
              $type .= "- tandem repeats"; 
              $countR->{tan}++;
            } elsif( $p2 < $o + $prtR[$l+2] ){
              $type .= "overlapping repeats"; 
              $countR->{ovr}++;
            } else {
              $type .= "dispersed repeats"; 
              $countR->{dis}++;
            }
            print MAP join("\t", @prtR[$l..$l+4], $type, $geneR),"\n";
            $type = "";
            $geneR = "";
          }
          @prtR = ();
          @tmpR = ();
          %seenR = ();

          $seenR{$pos} = 0;
          push @prtR, ($source, $ori."-".($ori+$kmer-1), $kmer, $r, ($pos+$kmer-1)."-".$pos);
          push @tmpR, ($source, $ori, $kmer, $r, $pos);
        }
      }
    }
  }
}
print "DONE\n";
close MAP;
close DIST;
$tree->clear();

if( $input_path ne ""){
  open HIST, ">$histfile";
  print HIST join("\t", "Dispersed", "Overlapping", "Tandem/Palindrome"),"\n";
  $ov = 0; $dis = 0; $ta = 0;
  $ov = $countF->{"ovr"};
  $ta = $countF->{"tan"};
  $dis = $countF->{"dis"};
  print HIST join("\t", "Direct", $dis, $ov, $ta, 0), "\n";
  $ov = 0; $dis = 0; $pal = 0;
  $ov = $countR->{"ovr"};
  $pal = $countR->{"pal"};
  $dis = $countR->{"dis"};
  $ta = $countR->{"tan"};
  print HIST join("\t", "Inverted", $dis, $ov, $pal, $ta), "\n";
  close HIST;
}

#################
## SUBROUTINES ##
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

