#!/usr/bin/perl

=head1 NAME

build_consensus_from_variants.perl

=head1 VERSION

0.1

=head1 DESCRIPTION

Using a range it returns a subseqeunce of a FASTA sequence

=head1 INPUT

FASTA file

=head1 OUTPUT

FASTA file

=head1 OTHER REQUIREMENTS

NONE

=head1 EXMAPLE USAGE

build_consensus_from_variants -i myseq.mpileup -r ref.fasta -o sample.fasta -l 20 -u 80 -c 20

=head1 AUTHOR

Aengus Stewart

=cut



use strict;
no warnings 'qw';
use Getopt::Std;
use FindBin qw( $Bin );

use vars qw( $opt_i $opt_r $opt_o $opt_l $opt_u $opt_h $opt_d $opt_c );
our ( $VERSION, @FIELD_NAMES );

$VERSION = 0.1;
@FIELD_NAMES = qw( CHR POS REF_BASE DEPTH SEQ BASE_Q MAP_Q );

my $correct_usage = "USAGE: build_consensus_from_variants [ flags ]\n\n" .
                    "  FLAGS\n" .
                    "    Required:\n" .
                    "      -i <sample MPILEUP>\n" .
                    "      -r <reference FASTA>\n" .
                    "      -l <num> lower threshold for ambiguity code\n" .
                    "      -u <num> upper threshold for ambiguity code\n" .
                    "      -c <num> depth coverage cutoff below which reference will be reported\n" .
                    "      -o <consensus FASTA>\n" .
                    "    Optional:\n" .
                    "      -d [1|2] debug level\n" .
                    "      -h help\n" .
                    "    It will incorporate an INDEL at a position if it occurs at >= 50%.\n" .
                    "    It will report in the log file positions where the depth is less then 10.\n" .
                    "    It will report in the log file positions where an INDEL occurs at >= 10%.\n";

my $help = "This will call a consensus from an mpileup file.\n";
my %amiguity_codes = ( "A" => "A",
                       "C" => "C",
                       "G" => "G",
                       "T" => "T",
                       "AC" => "M", 	
                       "AG" => "R",
                       "AT" => "W",
                       "CG" => "S",
                       "CT" => "Y",
                       "GT" => "K",
                       "ACG" => "V",
                       "ACT" => "H",
                       "AGT" => "D",
                       "CGT" => "B",
                       "GATC" => "N" );
my $record;
my $chr;
my $position;
my $ref_base;
my $depth;
my $seq;
my $base_q;
my $map_q;
my $R_seqs_reference;
my $R_seqs_chr;
my %seqs_out;
my @seq_tmp;
my $base;
my $base_count;
my $base_next;
my $base_percentage;
my $ambiguity = "";
my $length_indel;
my $length_deletion;
my $insertion;
my $ins_percentage;
my $del_percentage;
my $char;
my $i;
my $j;
my $j_offset;
my $k;
my $debug;

# Get command line args
getopts("i:o:r:u:l:d:c:h");

if ( ! $opt_i )
{
  print STDERR "PANTS: No input mpileup given\n\n$correct_usage\n\n\n";
  die;
}
else
{
  check_infile( *F_IN_MPILE, "${opt_i}" );
  check_outfile( *F_OUT_LOG, "${opt_i}.log" );
}

if ( ! $opt_o )
{
  print STDERR "PANTS: No output given\n\n$correct_usage\n\n\n";
  die;
}
else
{
  check_outfile( *F_OUT_FASTA, "${opt_o}" );
}

if ( ! $opt_r )
{
  print STDERR "PANTS: No input reference given\n\n$correct_usage\n\n\n";
  die;
}
else
{
  $R_seqs_reference = read_multiple_fasta_file( $opt_r );
}

if ( $opt_c !~ /\d+/ )
{
  print STDERR "PANTS: A depth coverage has to be given\n\n$correct_usage\n\n\n";
  die;
}
if ( $opt_l !~ /\d+/ )
{
  print STDERR "PANTS: A lower threshold has to be given for reporting an IUPAC ambuguity char\n\n$correct_usage\n\n\n";
  die;
}
if ( $opt_u !~ /\d+/ )
{
  print STDERR "PANTS: An upper threshold has to be given for reporting an IUPAC ambuguity char\n\n$correct_usage\n\n\n";
  die;
}

$debug = $opt_d if ( $opt_d );

#for my $key ( keys( %{ $R_reference_seqs } ) )
#{
#  print STDERR "$key\t" . $R_reference_seqs->{ $key} . "\n";
#  for my $key2 ( keys( %{ $R_reference_seqs->{ $key } } ) )
#  {
#    print STDERR "  $key2\n";
#    if ( $key2 =~ /seq/i )
#    {
#      print STDERR "    |" . $R_reference_seqs->{ $key }{ $key2 } . "|\n";
#    }
#  }
#}


while( defined( $record = <F_IN_MPILE> ) )
{
  # If we have a DELETION at > 50% need to skip these bases.
  if ( $length_deletion > 0 )
  {
    $length_deletion--;
    next;
  }

  chomp( $record );
  ( $chr, $position, $ref_base, $depth, $seq, $base_q, $map_q ) = split( /\s+/, $record );
  print STDERR "$chr $position\n" if ( $debug >= 1 );

  # If depth is low set to reference base
  goto GOTO_1 if ( $depth < $opt_c );

  # Parse the SEQ
  # Remove all START chars and the QUAL SCORE char following
  $seq =~ s/\^.//g;
  # Remove all STOP chars
  $seq =~ s/\$//g;
  #print STDERR "$seq\n";


  # Count the reference bases
  $base_count = 0;
  foreach $char ( qw(. ,) )
  {
    $base_count += length( $seq =~ s/[^\Q$char\E]//rg );
    print STDERR "$base_count\n" if ( $debug >= 2 );
    $seq =~ s/\Q$char\E//g;
  }
  $R_seqs_chr->{ $chr }[$position]{"\U$ref_base"} += $base_count;
  print STDERR "$seq\n" if ( $debug >= 2 );

  # We need to deal with the indels one char at a time.
  @seq_tmp = ( split( //, $seq ) );
  print STDERR "@seq_tmp\n" if ( $debug >= 2 );

  $i=0;
  while ( $i <= $#seq_tmp )
  {
    print STDERR "|$seq_tmp[$i]|\n" if ( $debug >= 2 );
    # Sort out DELETIONS - INCREASE DEL count and DECREASE BASE count
    if ( $seq_tmp[$i] =~ /[\+-]/ )
    {
      # Get the length of the indel
      $length_indel = "";
      $j_offset = 0;
      for $j ( ($i+1 ) .. $#seq_tmp )
      {
        $j_offset++;
        $length_indel .= $seq_tmp[$j] if ( $seq_tmp[$j] =~ /\d/ );
        if ( $seq_tmp[$j] !~ /\d/ )
        {
          $j_offset--;
	  last;
        }
      }
      print STDERR "length_indel: $length_indel\n" if ( $debug >= 1 );

      if ( $seq_tmp[$i] =~ /-/ )
      {
        $R_seqs_chr->{ $chr }[$position]{D}++;
        for $k ( 1 .. $length_indel )
        {
          # The next base will be the next char after the length identifier
          $base_next = $seq_tmp[($i+$j_offset+$k)];
          #print STDERR "del i=$i j=$j_offset k=$k position=$position |$base_next|\n";
          $R_seqs_chr->{ $chr }[($position+$k)]{$base_next}--;
          $R_seqs_chr->{ $chr }[$position]{DELETION_SEQ} .= $seq_tmp[($i+$j+$k)];
        }
      }
      elsif ( $seq_tmp[$i] =~ /\+/ )
      {
        $R_seqs_chr->{ $chr }[$position]{I}++;
        for $k ( 1 .. $length_indel )
        {
          $base_next = $seq_tmp[($i+$j_offset+$k)];
          #print STDERR "ins i=$i j=$j_offset k=$k position=$position |$base_next|\n";
          $R_seqs_chr->{ $chr }[$position]{INSERTION_SEQ} .= $seq_tmp[($i+$j+$k)];
        }
      }
      $i += $j_offset+$length_indel+1;
    }
    # This is ALTERNATIVE BASES
    else
    {
      $base = "\U$seq_tmp[$i]";
      print STDERR "ALT BASE = |$base|\n" if ( $debug >= 2 );
      $R_seqs_chr->{ $chr }[$position]{$base}++;
      $i++;
    }
  } # END while $i <= $#seq_tmp

  ######################################################################
  #
  # This has ended the parsing of the mpileup bases
  #
  ######################################################################

  GOTO_1:
  $ambiguity = "";
  print F_OUT_LOG "$chr " . ( $position ) . " Depth is low $depth\n" if ( $depth <= 10 );

  if ( $depth == 0 )
  {
    $del_percentage = $R_seqs_chr->{ $chr }[$position]{ D } / 1;
    $ins_percentage = $R_seqs_chr->{ $chr }[$position]{ I } / 1;
  }
  else
  {
    $del_percentage = $R_seqs_chr->{ $chr }[$position]{ D } / $depth;
    $ins_percentage = $R_seqs_chr->{ $chr }[$position]{ I } / $depth;
  }

  # Set BASE for current position
  if ( $depth <= $opt_c )
  {
    #$seqs_out{ $chr } .= $ref_base;
    if ( $depth == 0 )
    {
      $seqs_out{ $chr } .= "-";
    }
    else
    {
      $seqs_out{ $chr } .= "N";
    }
  }
  else
  {
    foreach $base ( qw( A C G T ) )
    {
      $base_percentage = $R_seqs_chr->{ $chr }[$position]{ $base } / $depth;
      print STDERR "$base " . $R_seqs_chr->{ $chr }[$position]{ $base } . " $base_percentage $depth\n" if ( $debug >= 2 );
      {
        # Homozygous
        if ( $base_percentage >= $opt_u )
        {
          $seqs_out{ $chr } .= $base;
          last;
        }
        # Heterozygous
        elsif ( $base_percentage >= $opt_l && $base_percentage <= $opt_u )
        {
          $ambiguity .= $base;
        }
      }
    } # END foreach
  } # END if ( $depth <= $opt_c )

    if ( $ambiguity ne "" )
    {
      $seqs_out{ $chr } .= $amiguity_codes{ $ambiguity };
    }
  
    # Deal with indels AFTER the current base is set as it refers to the NEXT base(s)
    if ( $del_percentage >= 0.5 )
    {
      print F_OUT_LOG "$chr " . ( $position + 1 ) . " deleted as deletion " . ( $del_percentage * 100 ) . "%\n";
      $seqs_out{ $chr } .= $R_seqs_chr->{ $chr }[$position]{ INSERTION_SEQ };
      $length_deletion = $length_indel;
      next;
    }
    elsif ( $del_percentage >= 0.1 )
    {
      print F_OUT_LOG "$chr " . ( $position + 1 ) . " deletion  " . ( $del_percentage * 100 ) . "%\n";
    }
  
    if ( $ins_percentage >= 0.5 )
    {
      print F_OUT_LOG "$chr " . ( $position + 1 ) . " insertion as insertion " . ( $ins_percentage * 100 ) . "%\n";
      $seqs_out{ $chr } .= $R_seqs_chr->{ $chr }[$position]{ INSERTION_SEQ };
    }
    elsif ( $ins_percentage >= 0.1)
    {
      print F_OUT_LOG "$chr " . ( $position + 1 ) . " insertion " . ( $ins_percentage * 100 ) . "%\n";
    }
      
    # RESET THESE VARS!!!!!!!
    $chr = "";
    $position = "";
    $ref_base = "";
    $depth = "";
    $seq = "";
    $base_q = "";
    $map_q = "";
  } # END while

foreach $chr ( keys( %seqs_out ) )
{
  print F_OUT_FASTA ">${chr}\n";
  print_fixed_lines( *F_OUT_FASTA, 50, $seqs_out{ $chr } );
}



sub check_infile
{

# Check if we can open the file 

#  print "$_[1]";
  open( $_[0], "<$_[1]" )
    or die "PANTS: Could not open input file $_[1]\n";
} # END check_infile


sub check_outfile
{

# Check if we can open the file

#  print "$_[1]";
  open( $_[0], ">$_[1]" )
    or die "PANTS: Could not open output file $_[1]\n";
} # END check_outfile

sub print_fixed_lines
{
  my ( $out_fh, $line_len, $stuff ) = @_;
  my $pos = 0;
 
  while ( $pos < length( $stuff ) )
  {
    print $out_fh substr( $stuff, $pos, ( $line_len ) ), "\n";
    $pos += ( $line_len );
  }
} # END sub print_fixed_lines

sub read_multiple_fasta_file
{
  my $fasta_file = shift;
  my $record;
  my $name;
  my $def;
  my %seqs;
  my $seq_dup_found;

  check_infile( *DATA_IN, "${fasta_file}" );

  while( defined( $record = <DATA_IN> ) )
  {
    chomp( $record );
    if ( $record =~ /^>/ )
    {
      $seq_dup_found = 0;
      ( $name, $def ) = $record =~ /^>\s*(\S+)\s*(.*)/;

      # Is there a valid sequence name
      if ( ! $name )
      {
        print STDERR "\nPANTS: No name for sequence in file $fasta_file\n";
        die;
      }

      # Have we already seen this sequence
      if ( ! defined( $seqs{ $name } ) )
      {
        $seqs{ $name }{ DEF } = $def;
        $seqs{ $name }{ SEQ } = undef;
      }
      else
      {
        print STDERR "WARNING: Duplicate sequence found will be ignored! - $name\n";
        $seq_dup_found = 1;
      }
    } # END if $record
    elsif ( ! $seq_dup_found )
    {
      # This is raw sequence
      $seqs{ $name }{ SEQ } .= $record;
    }
  }

  return( \%seqs );

} # END read_multiple_fasta_file



exit;

################################################################################
#
# Created                                                             03-12-2018
# Changed > $opt_c to N rather than REF BASE                          13-08-2019
# Changed $depth == 0 to report '-' rather than REF BASE              13-08-2019
