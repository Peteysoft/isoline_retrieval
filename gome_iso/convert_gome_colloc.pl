#!/usr/bin/perl

#for ($i=0; $i<scalar @ARGV; $i++) {
$dummy=0;
$nooverwrite=0;
for ($i=0; (substr($ARGV[$i], 0, 1) eq "-"); $i++) {
  #print substr($ARGV[$i], 0, 1), "\n"; 
  #if (substr($ARGV[$i], 0, 1) eq "-") {
    $opt=substr($ARGV[$i], 1);
    if ($opt eq "u") {
      $dummy=1;
    } elsif ($opt eq "n") {
      $nooverwrite=1;
    } else {
      print "Unknown command option:", $ARGV[$i], "\n";
    }
  #} else {
  #  break;
  #}
}

$argoffset=$i;

#print $argoffset;

$inputfilename=$ARGV[$argoffset];
$outpath=$ARGV[$argoffset+1];

#print "Opening: ", $inputfilename, "\n";
open INFILE, $inputfilename;

$fpref="coords";

$tempfile="temp";
$inbase="/misc/raider/GOME_LVL10/REPROC/";
$convert1="/misc/raider/GOME_LVL10/gdp01_ex_lx -b nynyyynnnn ";
$convert2="/home/home/pmills/my_software/isoret/convert_gome -dp ";

while (!eof(INFILE)) {
  $junk=<INFILE>;
  #print $junk;
  if (eof(INFILE)) {break;}

  $filename=<INFILE>;
  chop $filename;

  $junk=<INFILE>;
  #print $junk;

  #parse the filename:
  @path=split('\/', $filename);

  $np=scalar @path;

  $fbase=substr($path[$np-1], length($fpref));
  @suffixes=split('\.', $fbase);
  #print $fbase, @suffixes, "\n";
  $nsuf=scalar @suffixes;
  $fbase2=join('.', @suffixes[0..$nsuf-2]);

  if (substr($path[$np-1], 0, length($fpref)) ne $fpref) {
    print "Prefix of input file not correct:", $fpref, ", ", $path[np-1], "\n";
    exit(1);
  }

  $year=substr($fbase, 0, 4);
  $month=substr($fbase, 4, 2);
  $day=substr($fbase, 6, 2);

  $infile=$inbase.$year."/".$month."/".$day."/".$fbase2;
  #print @path[0..$np-2], "\n";
  $outfile=$outpath.join("/", @path[0..$np-2])."/".$fbase;

  #  $command1=$convert1.$infile." ".$tempfile;
  #  print $command1, "\n";
  if (!($nooverwrite && -f $outfile)) {
    $command1=$convert1.$infile." ".$tempfile;
    print $command1, "\n";
    if (!$dummy) {system $command1;}
    $command2=$convert2.$tempfile.".el1 ".$outfile;
    print $command2, "\n";
    if (!$dummy) {system $command2;}
  }
}

close(INFILE);
