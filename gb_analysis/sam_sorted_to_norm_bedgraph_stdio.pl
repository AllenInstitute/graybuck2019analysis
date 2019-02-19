#!/usr/bin/perl
# ---------------
# sam_sorted_to_norm_bedgraph_stdio.pl - converts a sorted, paired-end SAM file to a bedGraph file of fragment
# overlap counts normalized by the user-supplied fragment count argument. 
# Requires the read length of the raw reads used to generate the SAM file, as well as
# the minimum number of overlapping fragments you want displayed. The latter can greatly reduce the file
# size by leaving out orphaned fragments while retaining regions that are significantly enriched. Use 0
# for min_fragments to keep all reads. For this version, these values are set in the script, but could
# easily be changed for use as arguments.
# ---------------
# usage: samtools view input_file.bam | sam_sorted_to_norm_bedgraph_stdio.pl fragment_count > output_file.bedgraph
# ---------------
# SAM files must be from paired-end data, and can be generated from raw reads using the read alignment program Bowtie (http://bowtie-bio.sourceforge.net/).
# SAM files must be sorted before this script is used. SAM files can be sorted using IGVTools (http://www.broadinstitute.org/igv/igvtools).
# For more information about the SAM file format, see the SAMtools page (http://samtools.sourceforge.net/).
# For a description of bedGraph format, see the UCSC Genome Browser description of the format (http://genome.ucsc.edu/goldenPath/help/bedgraph.html).
# ---------------


$read_length = 50 - 1;
$min_frags = 0;

$frags = @ARGV[0];
$m_frags = $frags / 1000000;

$min_val = 1 / $m_frags;

print "track type=bedGraph\n";

$last_end = 0;
@cur_peak = ();
$cur_chr = "chrZ";

#subroutine for calculating the overlaps given a cur_peak array
sub calc_overlaps {
	#only process if there are at least the minimum number of fragments in the peak
	if(@cur_peak == 2) {
		
		($out_start, $val) = split(/\t/,@cur_peak[0]);
		($out_end, $val) = split(/\t/,@cur_peak[1]);

		print $cur_chr . "\t" . $out_start . "\t" . $out_end . "\t" . $min_val . "\n";
		
	} elsif( (@cur_peak > 2) && (@cur_peak >= $min_frags) ) {
		
		@sorted_peak = sort { $a <=> $b } @cur_peak;
		@cur_peak = ();
		
		($cur_start, $cur_val) = split(/\t/,@sorted_peak[0]);
					
		for ($i = 1; $i < @sorted_peak; $i++) {
			
			($pos, $val) = split(/\t/,@sorted_peak[$i]);
			
			if($val != 0) {
				
				$out_val = $cur_val / $m_frags;

				print $cur_chr . "\t" . $cur_start . "\t" . $pos . "\t" . $out_val . "\n";
				
				$cur_start = $pos;
				$cur_val += $val; 

			}
			
		}
		
		@sorted_peak = ();
		
	}	
}

while(<STDIN>) {

	$line = $_;
	chomp($line);
	@line = split(/\t/,$line);
	
	#Make sure this isn't a header line or from chromosome M
	if ( @line[0] eq "\@HD" || @line[0] eq "\@SQ" || @line[0] eq "\@PG" || @line[2] eq "chrM") {
	} elsif ( @line[1] == 99 || @line[1] == 163) { #if not, only process lines with flags 99 or 163.
	
		#get information about the fragment position
		$start = @line[3];
		#have to add the read length to the 3' end, because SAM format uses only the 5'-most position
		$end = @line[7] + $read_length;
		
		#Check to see if this is a new chromosome
		if ( @line[2] ne $cur_chr ) {
			#if so, output the final peak from the end of the last chromosome
			if($cur_chr eq "chrZ") { } else {
			calc_overlaps();
			}
			
			#and start processing the new chromosome.
			$cur_chr = @line[2];
			
			#reset the end point position and current peak array
			$last_end = 0;
			@cur_peak = ();
		
		} else {			
			#if this fragment is outside the last peak, output the last peak and start storing
			#fragments for a new peak.
			if($start > $last_end) {
			
				calc_overlaps();
			
				@cur_peak = ();
			}
			
		}
		
		#compare this read to the current peak array.
		
		$start_assigned = 0;
		$end_assigned = 0;
		
		if(@cur_peak > 0) {
			for($i = 0; $i < @cur_peak; $i++) {
				($pos,$val) = split(/\t/,@cur_peak[$i]);
				if ($start == $pos) {
					$val++;
					$cur_peak[$i]=$pos . "\t" . $val;
					$start_assigned++;
				} elsif ($end == $pos) {
					$val--;
					$cur_peak[$i]=$pos . "\t" . $val;
					$end_assigned++;
				}
			}
			
			if ($start_assigned == 0) {
				$new_start = $start . "\t1";
				push(@cur_peak, $new_start);
			}
			
			if($end_assigned == 0) {
				$new_end = $end . "\t-1";
				push(@cur_peak, $new_end);
			}
			
		} else {
			$new_start = $start . "\t1";
			push(@cur_peak, $new_start);
			$new_end = $end . "\t-1";
			push(@cur_peak, $new_end);
		}
					
		if($end > $last_end) {
			$last_end = $end;
		}
	}
}

#output the last peak

calc_overlaps();

close(SAM);

