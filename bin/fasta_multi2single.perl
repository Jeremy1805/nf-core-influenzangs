#!/usr/bin/perl

=head1 NAME

fasta_multi2single.perl

=head1 VERSION

0.1

=head1 DESCRIPTION

Split a multiple fasta file into single sequences

=head1 INPUT

FLAGS:
    B<-i> FASTA input file
    B<-o> FASTA basename for output files
    B<-d> debug level
    B<-h> (optional) help

=head1 OUTPUT

fasta files

=head1 OTHER REQUIREMENTS

NONE

=head1 EXAMPLE USAGE

fasta_multi2single.perl -i mybig.fasta -o mylittle

=head1 AUTHOR

Aengus Stewart CRUK

=cut

# START OF CODE

use strict;
use Fcntl qw(:seek);
use Getopt::Std;
use Time::HiRes qw(gettimeofday);
use FindBin qw( $Bin );

use vars qw( $opt_d $opt_h $opt_i $opt_o $VERSION );

$VERSION = 0.1;

my $count=0;
my $record;
my $file_count;
my $seq_name;
my $seq_def;
my $seq_seq;
my $seq_length;
my $ofilename;
my $found;
my $debug;
my $ofile_ext = "fa";

my $correct_usage = "USAGE: fasta_multi2single.perl [ flags ]\n\n" .
                    "  FLAGS\n" .
                    "    Required:\n" .
                    "      -i <fasta file>\n" .
                    "    Optional:\n" .
                    "      -o <string> basename for output files\n" .
                    "      -d debug <number>\n" .
                    "      -h\n";

my $help = "This program will do some stuff and then some other stuff\n";

# Get command line args
getopts("d:f:i:j:o:h");

# Set debug
$debug = $opt_d if ( $opt_d );

# Print out help
if ( $opt_h )
{
  print "\n$help\n\n$correct_usage\n\n\n";
  exit;
}

# Open the input file
die "No FASTA input file given!\n\n$correct_usage\n\n\n" if ( ! $opt_i );

print STDERR "Input          : $opt_i\n";
check_infile( *F_IN, "${opt_i}" );

# Check for output NAME
die "No Output Basename given! n\n$correct_usage\n\n\n" if ( ! $opt_o );



while ( defined( $record=<F_IN> ) )
{
  if ( $record =~ /^>/ )
  {
    $found = 1;
    $count++;
    $file_count = 1;
    $seq_length = length( $seq_seq );
    print_fixed_lines( \*SEQ_FILE, 50, $seq_seq );
    print "${opt_o}_${seq_name}\t$seq_length\n" if ( $seq_name ); 
    close( SEQ_FILE );
    $seq_seq = "";
    ( $seq_name ) = $record =~ /^>(.+?)[\s;]:*/;
    ( $seq_def ) = $record =~ /^>.+?([\s;]:*)/;
    $seq_name =~ s/[\, \`\'\"\(\)]/_/g;
    if ( $opt_o )
    {
      $ofilename = ${opt_o} . "." . $seq_name . "." . $ofile_ext;
    }
    else
    {
      $ofilename = $seq_name . "." . $ofile_ext;
    }
    $record = ">${opt_o}_${seq_name}${seq_def}";
    until ( ! $found )
    {
      if ( -e "${ofilename}" )
      {
        print "${ofilename}\n";
        if ( $opt_o )
        {
          $ofilename = ${opt_o} . "." . $seq_name . "." . $file_count . "." . $ofile_ext;
        }
        else
        {
          $ofilename = ${opt_o} . "." . $seq_name . "_" . $file_count . "." . $ofile_ext;
        }

        $record = ">${opt_o}_${seq_name}_${file_count}${seq_def}";
        $file_count++;
      }
      else
      {
        $found = 0;
      }
    }
    open( SEQ_FILE, ">${ofilename}" )
      or die "PANTS: Couldnt open ${ofilename}\n";
    print SEQ_FILE "$record";
  }
  else
  {
    chomp( $record );
    $seq_seq .= $record;
    #print SEQ_FILE $record;
  }
}

print_fixed_lines( \*SEQ_FILE, 50, $seq_seq );
print STDERR "Sequences found: $count\n";

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

sub check_infile
{

# Check if we can open the file 

#  print "$_[1]";
  open( $_[0], "<$_[1]" )
    or die "PANTS: Could not open input file $_[1]\n";
} # END check_infile



