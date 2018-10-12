#!/usr/bin/perl

=pod

=head1 Filter_Fastq_On_Mean_Quality


Filter_Fastq_On_Mean_Quality - Filter FASTQ reads on mean quality

=head1 SYNOPSIS

Filter_Fastq_On_Mean_Quality.pl -a ascii(default 33) -m mean_quality_threshold(default 30) -l length(default 35) -f fastq_file -o output_file 

=cut

use strict;
use Carp qw (cluck confess croak);
use warnings;
#use Readonly;
use Getopt::Long;
use Pod::Usage;
use Fatal qw/:void open close/;

=pod

=head1 OPTIONS

#--- describes parameters given to the script
#+++ command line syntax
#--- requirement of the option and its parameter can be:
#--- required: name or nature inside <>
#--- optional: name or nature inside []
#--- alternative between 2 elements: elements separated by a |

=head2 Parameters

=over 4

=item B<[option_name]> ([option nature]): #+++

[option description]. #+++
Default: [option default value if one] #+++

=back
#--- Example:
#---
#--- Template.pl [-help | -man]
#---
#--- Template.pl [-debug [debug_level]] [-size <width> [height]]
#---
#--- =over 4
#---
#--- =item B<-help>:
#---
#--- Prints a brief help message and exits.
#---
#--- =item B<-man>:
#---
#--- Prints the manual page and exits.
#---
#--- =item B<-debug> (integer):
#---
#--- Executes the script in debug mode. If an integer value is specified, it will
#--- be the debug level. If "-debug" option was used without specifying a debug
#--- level, level 1 is assumed.
#--- Default: 0 (not in debug mode).
#---
#---=item B<-size> (positive_real) (positive_real):
#---
#--- Set the dimensions of the object that will be drawn. The first value is
#--- the width; the height is the second value if specified, otherwise it will
#--- assume height and width are equal to the first value.
#--- Default: width and height are set to 1.
#---
#---=back

=cut
unless (@ARGV)
{
  pod2usage(0);
  exit();
}
# options processing
my ($man, $help, $debug, $ascii, $mean, $length, $fastq, $output) = (0, 0, 0, 33, 30, 35);
# parse options and print usage if there is a syntax error.
#--- see doc at http://perldoc.perl.org/Getopt/Long.html
GetOptions("help|?"     => \$help,
           "man"        => \$man,
           "debug:i"    => \$debug,
           "ascii|a=i" => \$ascii, 
           "mean|m=i" => \$mean,
           "length|l=i" =>\$length,
           "fastq|f=s"     => \$fastq,
           "output|o=s"  => \$output 
) or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}


my $file_handle;
my $output_handle;
my $score;
if (open($file_handle, $fastq) && open($output_handle,">$output"))
{
  while(<$file_handle>)
  {
    my $read=$_ . <$file_handle> . <$file_handle>; 
    $_ = <$file_handle>;
    chomp;
    map{ $score += ord($_)-$ascii} split(""); 
    print {$output_handle} "$read$_\n" if ($score/length($_)>= $mean && length($_)>=$length) ;
    $score= 0;
  }
  close($file_handle);
  close($output_handle);
}
else
{
     confess "ERROR: failed to open file:\n$!\n";
}

exit(0);

=pod

=head1 DIAGNOSTICS

#--- Give and explain here every error message the the script may generate.
#--- Use: (W)=warning (optional), (D)=deprecation (optional),
#---      (S)=severe warning (mandatory), (F)=fatal error (trappable),
#---      (X)=very fatal error (non-trappable).
#--- example:
#---
#--- =over 4
#---
#--- =item *
#---
#--- Negative radius;
#--- (F) Can not draw a circle with a negative radius.
#---
#---=back #+++

=head1 AUTHORS

#+++ Valentin GUIGNON (CIRAD), valentin.guignon@cirad.fr

[author1_name (company), email]#+++

[author2_name (company), email]#+++

=head1 VERSION

Version [version.subversion.build] #+++

Date [DD/MM/YYY] #+++

=head1 SEE ALSO

#--- other documentation or objects related to this package
[perl(1), L<Geometry::Square>] #+++

=cut

