#This script obtained from Rahul Kulkarni University of Massachusetts, Boston,
#(https://www.umb.edu/academics/csm/faculty_staff/rahul_kulkarni)
#Some small modifications made to generalise parts of the code made by JPJH November 2016

#This program aims to find targets of the protein CsrA. The algorithm is described in the paper:
#http://nar.oxfordjournals.org/cgi/pmidlookup?view=long&pmid=24782516


use strict;
use FileHandle;

my $FASTA_START = -200;
my $FASTA_END = 30;
my $SD_REGION_LENGTH = 30;

#this is the main function that tracks possible csra targets#
#------------------------------------------------------------------------------#
sub csratargets
{
	my $orgName = $_[0];
	my $geneName;
	my $fastaSeq;
	my $subseq;                                        
	my $SDseq;                                         
	my $lengthOrf;										        									
	my $upstreamSeq;
	my $lengthUpstream;
	my @new_begin;
	my @new_end;
	my $new_positions;
    my $line;
	my $line2;
    my $totalTargets = 0;							
    my $temp1;
    my $temp2;
    my $temp;
	my $total_genes = 0;
    my $total_sd_present = 0;
	my $motifCount;
    my $motif;
    my $SD_found;
    my @begin_positions;
    my @end_positions;
	my $motif_position;
    my $i;			
    my $j;
	my $k;
	my $m;
	my $p;
	my $q;
	my $d1 = 10;
	my $d2 = 60;
	my $d3 = 10;
	my $tss_found;
	my $site_type;
	my $site_type_s;
	my $site_type_w;
	my $temp_position;
	my @temp_positions;
	my $temp_prev;
	my $temp_next;
	my @positions;
    my @motif_type;
    my $search_motif;
   	my $possible_targets;
	my @tss_details = ();
	my $tss_entries = 0;
    my $inFile;
	my $inFile2;
    my $outFile;
    my $outFile2;
    
    chomp $orgName;
    $inFile = $orgName.'_regions.txt';
	$inFile2 = $orgName.'_tss_details.txt';
	open(INFILE,$inFile) || die "sorry, system can't open $inFile.";
	if(-e $inFile2)
	{
		open(INFILE2,$inFile2) || die "Sorry, systen can't open $inFile2.";
		$line2 = <INFILE2>;
		while(!(eof(INFILE2)))
		{
			if($line2 =~/(\S+)\s+(\d+)/)
			{
				$tss_details[$tss_entries][0] = $1;
				$tss_details[$tss_entries][1] = $2;
				$tss_entries++;
				$line2 = <INFILE2>;
			}
			else
			{
				$line2 = <INFILE2>;
			}
		}
	}
	$outFile = "csra_targets_$orgName";
	open(OUTFILE1,">$outFile.txt");
    $outFile2 = "genenames_$orgName";
    open(OUTFILE2,">$outFile2.txt");
    print "Processing $orgName. Please wait...\n";
    print OUTFILE1 "Possible CsrA Targets:\n\n";
	print OUTFILE2 "gene\tNo. of primary binding sites\tNo. of secondary binding sites\n";
	$line = <INFILE>;
    while(!(eof(INFILE)))
    {
      if($line =~ /^\>(\S+)\s+(\S+);.*to \+{0,1}(\d+)/) # removed some of the conditions here
	# to allow input from other sources (JH Aug 2019)
    	{
			$total_genes++;
			$geneName = $1;
    		$lengthOrf = $3 + 1;
			$line = <INFILE>;
    		$fastaSeq = ();
			until (($line =~/^\>/) || eof(INFILE))
    		{
    			chomp($line);
    			$fastaSeq.= lc($line);									
    			$line = <INFILE>;
    		}
    		$lengthUpstream = length($fastaSeq) - $lengthOrf;
    		$upstreamSeq = substr($fastaSeq,0,$lengthUpstream+5);
    		$SDseq = substr($upstreamSeq,-$SD_REGION_LENGTH-5);
			$SD_found = 0;
    		while($SDseq =~ /(a.gga)/g)							
         	{
         	    $SD_found = 1;
         	    pos($SDseq)--;
			}
			
         	#if there is no aNgga found, find the last agga as the SD sequence#
         	if($SD_found != 1)
         	{
  	    
				while($SDseq =~ /(agga)/g)
         	    {
         	        $SD_found = 1;
         	        pos($SDseq)--;
        	    }
         	}
         	if($SD_found == 1)
         	{
				$total_sd_present++;
				$tss_found = 0;
				for($i=0;$i<$tss_entries;$i++)
				{
					if(($tss_details[$i][0] eq $geneName) && $tss_details[$i][1]<200)
					{
						$temp1 = $lengthUpstream-$tss_details[$i][1];
						$tss_found = 1;
						last;
					}
				}	
				if($tss_found ==0)
				{
					$temp1=$FASTA_START+$lengthUpstream;
				}
				$temp2 = $FASTA_END - $FASTA_START+1 ;
    			$subseq = substr($fastaSeq,$temp1,$temp2);
				$motifCount = 0; 
   				@begin_positions = (); @end_positions = ();
   				$motif_position = -1;
				$i=0;$temp_prev = 0;$temp_next = 0;$m=0;@temp_positions = 0;@motif_type = ();@positions = 0;
				while(($subseq =~/(.{0,4}(a|g)ga)/g))
				{
					$search_motif = $1;
					if(($search_motif =~ /(a(t|c|g|a){0,1}gga)/) ||  ($search_motif =~ /(ctgga)/) || ($search_motif =~ /(agaga)/) ||  ($search_motif =~ /(cggga)/) || ($search_motif =~ /(tggga)/))
					{
						$motif = $1;
						$motif_position = index($subseq,$motif,$motif_position+2);
						if($motif eq 'ctgga' || $motif eq 'agaga' || $motif eq 'cggga' || $motif eq 'tcgga')
						{
							$motif_type[$i] = 'w';
						}
						else
						{
							$motif_type[$i] = 's';
						}	
						$temp_positions[$i] = $motif_position;
						$i++;
						pos($subseq) = pos($subseq)- (length($search_motif) - index($search_motif,$motif)-length($motif) + 1);
					}	
					else
					{
						pos($subseq) = pos($subseq)-length($search_motif)+2;
					}
				}
				
				for($j=0;$j<$i;$j++)
				{
					if($motif_type[$j] eq 'w')
					{
						if($j==0)
						{
							$temp_prev = 0;
						}
						else
						{
							$temp_prev = $temp_positions[$j]-$temp_positions[$j-1];
						}
						if($j==$i-1)
						{
							$temp_next = 0;
						}
						else
						{
							$temp_next = $temp_positions[$j+1]-$temp_positions[$j];
						}
						if($temp_prev>$d3 && $temp_next>$d3)
						{
							next;
						}
						else
						{
							$positions[$m]=$temp_positions[$j];
							$m++;
						}
					}
					else
					{
						$positions[$m]=$temp_positions[$j];
						$m++;
					}
				}
				for($j=0;$j<$m-1;$j++)
				{
					$temp = $positions[$j+1] - $positions[$j];
				    if($temp>$d2)
					{
						$subseq = substr($subseq,$positions[$j+1]);
						$temp_position = $positions[$j+1];
						$positions[$j+1] = 0;
						for($k=$j+2;$k<$m;$k++)
						{
							$positions[$k] = $positions[$k]-$temp_position;
						}	
					}
				}
				@motif_type = ();$motif_position =0;
  				while(($subseq =~/(.{0,4}(a|g)ga)/g))
  				{
 					$search_motif = $1;
 					if(($search_motif =~ /(a(t|c|g|a){0,1}gga)/) ||  ($search_motif =~ /(ctgga)/) || ($search_motif =~ /(agaga)/) ||  ($search_motif =~ /(cggga)/) || ($search_motif =~ /(tggga)/))
 					{
						$motif = $1;
						if(($motif eq 'aagga') || ($motif eq 'atgga') || ($motif eq 'acgga') || ($motif eq 'aggga') || ($motif eq 'agga'))
   						{
  							$motif_type[$motifCount] = 's';
  						}
   						elsif(($motif eq 'ctgga') || ($motif eq 'cggga') || ($motif eq 'tggga') ||($motif eq 'agaga'))
   						{
  							$motif_type[$motifCount] = 'w';
   						}
   						$begin_positions[$motifCount] = index($subseq,$motif,$motif_position);
   						$end_positions[$motifCount] = $begin_positions[$motifCount]+length($motif)-1;
						$motif_position = $end_positions[$motifCount];
						$motifCount++;
    					pos($subseq) = pos($subseq)- (length($search_motif) - index($search_motif,$motif)-length($motif) + 1);
 					}
 					else
 					{
						pos($subseq) = pos($subseq)-length($search_motif)+2;
					}
  				}
  				if($motifCount >= 3)
  				{
 					$j=0;
 					@new_begin = 0;
 					@new_end = 0;
					$site_type_s = 0;
					$site_type_w = 0;
 					for($i=0;$i<$motifCount;$i++)
 					{
						$site_type = $motif_type[$i];
						if($i==0)
						{
						    $temp1 = 0;
   							$temp2 = $begin_positions[$i+1]-$end_positions[$i];
						}
						elsif($i==$motifCount-1)
						{	
						    $temp1 = $begin_positions[$i]-$end_positions[$i-1];
   							$temp2 = 0;
						}
						else
						{
						    $temp1 = $begin_positions[$i]-$end_positions[$i-1];
   							$temp2 = $begin_positions[$i+1]-$end_positions[$i];
						}	
						if($site_type eq 'w' && (($i==0 && $temp2<=$d3) || ($i==$motifCount-1 && $temp1<=$d3) ||     (($i>0 && $i<$motifCount-1) && ($temp1<=$d3 || $temp2 <=$d3))))
						{
						    $new_begin[$j] = $begin_positions[$i];
   							$new_end[$j] = $end_positions[$i];
   							$j++;
							$site_type_w++;
						}
						elsif($site_type eq 's')
						{
						    $new_begin[$j] = $begin_positions[$i];
   							$new_end[$j] = $end_positions[$i];
   							$j++;
							$site_type_s++;
						}
 					}		
 					$new_positions = $#new_begin+1;
					$possible_targets = 0;
					if($site_type_s >2 or ($site_type_s ==2 && $site_type_w >1))
					{
						for($i=0;$i<$new_positions-1;$i++)
						{
							for($j=$i+1;$j<$new_positions;$j++)
							{
								$temp = $new_begin[$j] - $new_end[$i];
								if($temp <0)
								{
									$temp = 0;
								}
								if($temp >=$d1 && $temp <=$d2)
								{
									$possible_targets++;
								}
							}
						}
					}
					if($possible_targets >1)
					{	
						print OUTFILE1 ">$geneName\n$subseq\n\n";
						print OUTFILE2 "$geneName\t$site_type_s\t$site_type_w\n";
						$totalTargets++;
						next;
					}	
				}
			}
        }
    	else
    	{
    		$line = <INFILE>;
    	}
	}			
	print "total genes:$total_genes\n";
    print "Genes in SD:$total_sd_present\n";
    print OUTFILE1 "\nTotal Entries:$totalTargets\n";
    print OUTFILE2 "\nTotal Entries:$totalTargets\n";
	print "\ntotal entries:$totalTargets\n\n";

    close INFILE;
	close INFILE2;
    close OUTFILE1;
    close OUTFILE2;

    print "Processing $orgName complete.\nPlease see $outFile and $outFile2 for output.\n\n";
}

#------------------------------------------------------------------------------#

#output function#

my $ORGFILE = new FileHandle ("organisms.txt") || die "can not open orgfile\n";
my $org;

while ($org = <$ORGFILE>)
{
    if($org =~ /(\S+)/)
    {
        $org = $1;
        &csratargets($org);
    }
}

close ($ORGFILE);

