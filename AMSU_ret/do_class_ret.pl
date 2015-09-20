#!/usr/bin/perl

$basepath="/freax/storage/home/pmills/isoline_ret2/";

$amsu_base="/talax/storage/noaa/";
$satid=15;

$level3_afile=$basepath."ret/a.l1c";
$level3_bfile=$basepath."ret/b.l1c";

$script="/freax/storage/home/pmills/software/atovs_tools/shell_scripts/zamsu2l1c.sh";

$year='2002';
$month='09';

$err=0;
for ($satid=15; $satid<=17; $satid++) {

for ($day=1; $day<=9; $day++) {
  $daystr=sprintf("%2.2d", $day);

  $amsua_path="noaa".$satid."_amsua_".$year."/".$month."/".$daystr."/";
  $amsub_path="noaa".$satid."_amsub_".$year."/".$month."/".$daystr."/";

  print $amsu_base.$amsua_path, "\n";

  open AF, "ls ".$amsu_base.$amsua_path."|";
  $j=0;
  while(<AF>) {
    chop;
    $afiles[$j]=$_;
    $j++
  }
  close AF;
  open BF, "ls ".$amsu_base.$amsub_path."|";
  $j=0;
  while(<BF>) {
    chop;
    $bfiles[$j]=$_;
    $j++
  }
  close BF;

  #print @afiles;

  $na=scalar(@afiles);
  print $na, "\n";
  for ($j=0; $j<$na; $j++) {
    $blen=length($afiles[$j]);
    $check_bfile=$afiles[$j];
    substr($check_bfile, 6, 1)="B";
    if ($bfiles[$j] ne $check_bfile) {
      printf, "A and B names do not agree: ", $afiles[$j], $bfiles[$j];
      stop;
    }
    $resultfile=$basepath."data/cret/".$year.$month.$daystr.
		"/CRET.".substr($afiles[$j], 9, 33).".r1c";

    #convert the level 1b files to level 1c:
    #shell_command =  'ssh marvin "' + script + " " + $
    #            amsu_base+amsua_path+afiles[$j] +" " + level3_afile + '"'
    $shell_command =  $script." ".$amsu_base.$amsua_path.$afiles[$j]." ".$level3_afile;
    print $shell_command, "\n";
    $err=system($shell_command);
    if ($err == -1 || $err >> 8 != 0) {
      print "Failed to convert file: ", $afiles[$j];
      exit;
      next;
    }
    #shell_command =  'ssh marvin "'.$script + " " + 
    #            amsu_base+amsub_path+bfiles[$j] +" ".$level3_bfile.'"';
    $shell_command =  $script." ".$amsu_base.$amsub_path.$bfiles[$j]." ".$level3_bfile;
    print $shell_command, "\n";
    $err=system($shell_command);
    if ($err == -1 || $err >> 8 != 0) {
      print "Failed to convert file: ", $bfiles[$j];
      exit;
      next;
    }
    
    $command=$basepath."src/retrieve_isoline2 ".$level3_afile." ".$level3_bfile." ".$resultfile;
    print $command, "\n";
    $err=system($command);
    if ($err == -1 || $err >> 8 != 0) {
      print "Retrieval for files ",$afiles[$j]," and ", $bfiles[$j], " failed.\n";
      exit;
    }

  }
}

}


