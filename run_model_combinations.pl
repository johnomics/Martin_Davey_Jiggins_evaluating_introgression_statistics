#!/usr/bin/env perl

# run_model_combinations.pl
# Run a set of combined models using shared_ancestry_simulator.R based on a CSV file of model parameter combinations

# Written for "Evaluating statistics for the identification of introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# John Davey:   jd626@cam.ac.uk
# Simon Martin: shm45@cam.ac.uk
# November-December 2013

use strict;
use warnings;
use Carp;
use English;
use Data::Dumper;
use Getopt::Long;


my $numwindows = 1000;
my $windowsize = 5000;
my $modelfile = "";
my $threads = 1;
my $popsize = 1000000;
my $simulator_path = "..";
my $options_okay = GetOptions(
	'numwindows=i' => \$numwindows,
	'windowsize=i' => \$windowsize,
	'modelfile=s' => \$modelfile,
	'threads=i'   => \$threads,
    'popsize=i'   => \$popsize,
    'simulator_path=s' => \$simulator_path,
);

croak "Please supply a CSV file containing model parameters with -m" if $modelfile eq "";

croak "Simulator path $simulator_path does not exist!\n" if !-d $simulator_path;

my $model = <<EOM;
seqgen:
    m: HKY
    l: $windowsize
ms:
    pops:
        - {name: 1, seqs: 8}
        - {name: 2, seqs: 8}
        - {name: 3, seqs: 8}
        - {name: 4, seqs: 8}
    opts:
        genperyear: 4
        n0: $popsize
        mu: 2.5e-9
        split:
            - {time: 3,   pops: [4,1]}
EOM

sub output_models {
    my ($model, $modeltype, $params) = @_;
    my $pops = $modeltype eq "Background" ? "2,1" : $modeltype eq "Alternate" ? "2,3" : "0,0";
    my $printpops = $pops;
    $printpops =~ s/,//g;
    foreach my $t123 (sort {$a<=>$b} keys %{$params}) {
        my $model_t123 = $model . "            - {time: $t123, pops: [3,1]}\n";
        foreach my $txx (keys %{$params->{$t123}}) {
            open my $model_file, ">", "$modeltype\_t123-$t123\_t$printpops-$txx.yml" or croak "Can't open model file: $OS_ERROR!\n";
            my $model_t123_txx = $model_t123 . "            - {time: $txx, pops: [$pops]}\n";
            print $model_file $model_t123_txx;
            close $model_file;
        }
    }
}

open my $combined, "<", $modelfile or croak "Can't open model file: $OS_ERROR!\n";

my $header = <$combined>;
chomp $header;
my @header = split /,/, $header;
my %background;
my %alternate;
my @combinations;
while (my $model = <$combined>) {
    chomp $model;
    my @f = split /,/, $model;
    $background{$f[1]}{$f[2]}++;
    $alternate{$f[3]}{$f[4]}++;
    my %comb;
	$comb{type}=$f[0];
    map {
        my ($type, $var) = split /_/, $header[$_];
        $comb{pars}{$type}{$var} = $f[$_]
    } 1..$#header;
    push @combinations, \%comb;
}

close $combined;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
$mon++; # originally in 0..11
$year+=1900;
my $outdir = "model_files_win$numwindows\_size$windowsize\_ne$popsize\_$year-$mon-$mday";
mkdir($outdir) if !-d $outdir;
chdir($outdir);

open my $log, ">", "$outdir.log" or croak "Can't open log file! $OS_ERROR\n";

output_models($model, "Background", \%background);
output_models($model, "Alternate", \%alternate);


foreach my $combination (@combinations) {
    my $command = "$simulator_path/shared_ancestry_simulator.R -w $numwindows -T $threads -c ";
    foreach my $model_type (keys %{$combination->{pars}}) {
        my $prop = $model_type eq "Background" ? 0.9 : 0.1;
        $command .= $model_type . "_t123-";
        $command .= $combination->{pars}{$model_type}{t123};
        $command .= "_";
        delete $combination->{pars}{$model_type}{t123};
        my $var2 = (keys %{$combination->{pars}{$model_type}})[0];
        $command .= "$var2-$combination->{pars}{$model_type}{$var2}.yml:$prop,";
    }
    chop $command;
    print $log "$command\t" . localtime() . "\n";
    my $simout = `$command`;
    croak "Couldn't run command: $command\n(Is shared_ancestry_simulator.R in your path?)\n" if !defined $simout;
    croak "Command failed: $command\n$simout\n" if ($simout ne "");
}