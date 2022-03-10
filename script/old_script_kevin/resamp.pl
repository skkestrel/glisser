#!/usr/bin/perl -w
use List::Util qw(first max min sum);
use POSIX;

######imporant notes about input file format:
# this file assumes that test particles have positive id numbers
# and planets have negative id numbers (swift default)
# and the default planet if none is specified is Neptune, id=-5

# the file is assumed to have the format:
# <id> <time> <semimajor axis> <eccentricity> <inclination(deg)> <long. of asc. node (deg)> <argument of perihelion (deg)> <mean anomaly (deg)>
#

# the libration characteristics are calculated over a set of fixed windows in time (could easily be altered
# to use sliding windows, but it doesn't matter for well-behaved libration

use Getopt::Long;
#arguments form the command line: minimum requirement of -if <input filename> -of <output filename> -tp <ntp> -p <p> -q <q>
my $opt_status = GetOptions("if=s" => \$infile, #follow file name
					 "of=s" => \$outfile, #follow file name
					 "tp=s" => \$ntp, #number of test particles
					 "p=s" => \$p, #res arg of the form: 
										# $p*(lambda_kbo) - $q*(lambda_planet) - $mr*(l-peri_kbo) - $nr*(node_kbo) 
										# - $rr*(l-peri_planet) - $sr*(node_planet);
					 "q=s" => \$q,
					
					#optional arguments below this line
					"w=s" => \$nwindows, #number of windows to divide up the simulation into, default is 20
					"m=s" => \$mr,
					"n=s" => \$nr,
					"r=s" => \$rr,
					"s=s" => \$sr,
					"pl=s" => \$pl); #id that defines planet (enter as positive number)

if(!$p or !$q or !$infile or !$outfile or !$ntp)
	{
	print "required arguments missing\n";
	die;
	}


#default to simplest eccentricity resonance
if(!$mr)
	{
	$mr = $p-$q;
	$nr=0;
	$rr=0;
	$sr=0;
	}
if(!$nr){$nr=0;}
if(!$rr){$rr=0;}
if(!$sr){$sr=0;}

#default to Neptune
if(!$pl){$pl=5;} 

#default to 20 windows
if(!$nwindows){$nwindows=20;}

$pi = 3.141592653589793;


$pl = -$pl;



for($i=0;$i<=$ntp;$i++)
	{
	$particlefound[$i] = 0;
	$atot[$i] = 0;
	$etot[$i] = 0;
	$itot[$i] = 0;
	$ptot[$i] = 0;
	$tsmax[$i] = 0;
	}

$s = -1;
open(DAT,"<$infile");
while(<DAT>)
	{ #start reading file
	$stuff = $_;
	chomp($stuff);
	$stuff = " ".$stuff;
	@data = split(/\s+/,$stuff);
	$node = $data[6]*$pi/180;
	$peri = $data[7]*$pi/180;
	$ma = $data[8]*$pi/180;
	$l = $node+$peri+$ma;
	if($data[1] == $pl)
		{
		$s = $s+1;
		$lambda_N = $l;
		$node_N = $node;
		$peri_N = $peri;
		}
	if($data[1] >= 0)
		{
		$j = $data[1];
		#calculate the resonant argument
		$tp = $p*$l - $q*$lambda_N - $mr*($peri+$node) - $nr*$node - $rr*($peri_N+$node_N) - $sr*$node_N;
		while ($tp > 2.0*$pi){$tp = ($tp-2.0*$pi);}					
		while ($tp < 0.0){$tp = ($tp+2.0*$pi);}
		$ptot[$j] = $ptot[$j]+$tp;
		$psi[$j][$s] = $tp;
		
		$a[$j][$s] = $data[3];
		$ec[$j][$s] = $data[4];
		$inc[$j][$s] = $data[5];

		
		$atot[$j] = $atot[$j] + $data[3];
		$etot[$j] = $etot[$j] + $data[4];
		$itot[$j] = $itot[$j] + $data[5];
		$tsmax[$j] = $s;
		$particlefound[$j] = 1;
		}
	} #stop reading file
$nlines = $s;
$asmax = $nlines;

#set up output file for the resonance amplitudes
open(AMPOUT,">$outfile");
print AMPOUT "# id phibar, dphi; ebar, de; ibar, di; abar, da; a0, e0, i0; fraction of windows librating; fraction of sim;  phimax-phimin\n";
print AMPOUT "# 1;   2,     3;	  4,  5;    6,	 7;    8,   9; 10, 11, 12;                  13;                14;                15       \n";

#determine the window size to use for the amplitude calculations:
$pwindow = floor($nlines/$nwindows);
$dpmax = 0;
$dpmin = 300;

for($n=0;$n<=$ntp;$n++)
	{ #start loop over test particles
	if ($particlefound[$n] == 0)
		{
		next;
		}
	$smax = $tsmax[$n];
	$s = 0;
	$abar = $atot[$n]/($smax+1.0);
	$pbar = $ptot[$n]/($smax+1.0);
	$pbar = $pbar*180/$pi;
	$ebar = $etot[$n]/($smax+1.0);
	$ibar = $itot[$n]/($smax+1.0);
	$dpavg = 0.0;
	$daavg = 0.0;
	$deavg = 0.0;
	$diavg = 0.0;
	$count = 0;
	$jmax = $smax-$pwindow;
	$lib_count = 0;
	for($j=0;$j<=$jmax;$j+=$pwindow)
		{
		$ik = $j+$pwindow;
		@temparray = ();
		@temparray = @{$psi[$n]}[$j..$ik];
		
		@sorted_p = ();
		@sorted_p = sort {$a <=> $b} @temparray;
		#take three measures of the libration amplitude and average them (averages over outliers)
		$p_window_amp = ($sorted_p[$pwindow]-$sorted_p[0]) + ($sorted_p[$pwindow-1]-$sorted_p[0]) +
								($sorted_p[$pwindow]-$sorted_p[1]);
		$p_window_amp = $p_window_amp/3.0;
				
		$p_window_amp = $p_window_amp*90/3.141592653589793;
		if($p_window_amp < 170){$lib_count = $lib_count+1;} #170 is an arbitrary choice that can break
																			 #down for some of the less well behaved resonances
		$dpavg = $dpavg + $p_window_amp;

		#calculate the changes in a,e,i over the window
		@temparray = ();
		@temparray = @{$a[$n]}[$j..$ik];
		@sorted_p = ();
		@sorted_p = sort {$a <=> $b} @temparray;
		$a_window_amp = ($sorted_p[$pwindow]-$sorted_p[0]) + ($sorted_p[$pwindow-1]-$sorted_p[0]) +
								($sorted_p[$pwindow]-$sorted_p[1]);
		$a_window_amp = $a_window_amp/3.0;
		$daavg = $daavg + $a_window_amp;

		@temparray = ();
		@temparray = @{$ec[$n]}[$j..$ik];
		@sorted_p = ();
		@sorted_p = sort {$a <=> $b} @temparray;
		$e_window_amp = ($sorted_p[$pwindow]-$sorted_p[0]) + ($sorted_p[$pwindow-1]-$sorted_p[0]) +
								($sorted_p[$pwindow]-$sorted_p[1]);
		$e_window_amp = $e_window_amp/3.0;
		$deavg = $deavg + $e_window_amp;

		@temparray = ();
		@temparray = @{$inc[$n]}[$j..$ik];
		@sorted_p = ();
		@sorted_p = sort {$a <=> $b} @temparray;
		$i_window_amp = ($sorted_p[$pwindow]-$sorted_p[0]) + ($sorted_p[$pwindow-1]-$sorted_p[0]) +
								($sorted_p[$pwindow]-$sorted_p[1]);
		$i_window_amp = ($i_window_amp/3.0)*180.0/$pi;
		$diavg = $diavg + $i_window_amp;

		$count = $count + 1;
		}
	$dp = $dpavg/$count;
	$da = $daavg/$count;
	$de = $deavg/$count;
	$di = $diavg/$count;

	$jkl = 0;
	@temparray = ();
	@temparray = @{$psi[$n]}[$jkl..$smax]; @sorted_p = ();
	@sorted_p = sort {$a <=> $b} @temparray;
	$tp = ($sorted_p[$smax]-$sorted_p[0])*90/$pi; 
	if($dp > $dpmax){$dpmax = $dp;}
	if($dp < $dpmin){$dpmin = $dp;}

	#output the results along with a fractional resonance libration stat			
	
	$temp = $lib_count/$count;
	$fint = $smax/$asmax;
	print AMPOUT "$n $pbar $dp $ebar $de $ibar $di $abar $da $a[$n][0] $ec[$n][0] $inc[$n][0] $temp $fint $tp\n";
	} #end loop over test particles

close(AMPOUT);

