#!/usr/bin/env perl

# Author: Murray Cox

# Programme: mstree v.2.0
# A programme to parse tree info from ms data
# Institution: University of Arizona and Santa Fe Institute
# Date: May 2007
# Developed on perl v5.8.6 built for darwin-thread-multi-2level

# Checked that code still runs
# Date: December 2018
# Confirmed on perl v5.18.2 built for darwin-thread-multi-2level

# This is a release version of mstree.pl
# Print statements for non-essential statistics have been excised
# However, code to calculate non-essential statistics is still executed
# Release code not optimized for runtime speed

# example: ms 10 1 -t 5 -I 2 6 4 1 -T | mstree -g1 6 -g2 4

use strict;
use warnings;
use Getopt::Long;

# flag values
my( $group_one, $group_two );
GetOptions(
	"g1=s" => \$group_one,	# group one sum
	"g2=s" => \$group_two	# group two sum
); 

# usage information
my $usage = "Error: correct usage is: mstree -g1  number of individuals in group one
                                -g2  number of individuals in group two
                                [group 1 must be the selected group]\n";

die $usage if ( !$group_one || !$group_two || $group_one < 1 || $group_two < 1 );

# tmrca in scalar x Ne generations
print STDOUT "pmc\ttmrca\n";

# global variables
my $past_first_segsites_line = 0;
my( $ms_size, $segsites, @sequences, $r_prop, $m_prop, @root_n_persistent, @main_n_persistent );

MAIN: while ( <STDIN> ) {
	
	# count variables
	my $lb_c = 0;
	my $rb_c = 0;
	my $n_count = 0;
	my $l_count = 0;
	
	# storage arrays
	my( @root_n, @root_l, @main_n, @main_l );
	my( @bracket_order );
	
	# read in line
	my $line = $_;
	chomp($line);
		
	# check for ms line
	if ( $line =~ /ms (\d+)/ ) {
		$ms_size = $1;
		if ( $ms_size != ( $group_one + $group_two ) ) {
			die "error: input flags do not sum to $ms_size\n";
		}
	}
	
	# count number of segregating sites
	if ( $line =~ /segsites: (\d+)/ ) {
		$past_first_segsites_line = 1;
		$segsites = $1;
	}
	
	# extract to sequence matrix if line starts with 0 or 1
	if ( $past_first_segsites_line == 1 && $line =~ m/^(0|1)/i ) {
		chomp($line);
		push(@sequences, $line);
	}

	# sequence analysis
	if ( scalar(@sequences) == $ms_size ) {
		
		# determine minimum proportion clade and assign hashes
		my( %min_clade, %max_clade );
		if( $r_prop < $m_prop ){
			
			foreach my $root1 (@root_n_persistent){
				$min_clade{$sequences[$root1 - 1]} = "MIN";
			}
			foreach my $main1 (@main_n_persistent){
				$max_clade{$sequences[$main1 - 1]} = "MAX";
			}
		} elsif( $m_prop < $r_prop ){
			
			foreach my $main1 (@main_n_persistent){
				$min_clade{$sequences[$main1 - 1]} = "MIN";
			}
			foreach my $root1 (@root_n_persistent){
				$max_clade{$sequences[$root1 - 1]} = "MAX";
			}
		} 
		
		# merge hashes
		my %haplotypes = (%min_clade, %max_clade);
		
		# extract number of unique lineages in pmc clade
		my @unique = grep {$_ =~ m/^MIN/} values %haplotypes;
		my $pmc_haplotypes = scalar(@unique);
		
		# clean up variables
		@sequences = ();
	}
		
	# ignore non-tree lines
	next MAIN if $line !~ /;/;

	# collect characters
	my @chars = split(//, $line);
	
	# remove first and last two characters
	shift @chars;
	pop @chars; pop @chars;
	
	# count parentheses (if first state is a parenthesis)
	if ( $chars[0] eq '(' ) {
		BRACKETS: foreach my $entry ( @chars ) {
			
			$lb_c++ if $entry eq '(';
			$rb_c++ if $entry eq ')';
			
			last BRACKETS if $lb_c == $rb_c;
		}
	}
	
	# record order of parentheses (for all states)
	my $lp_c = 0;
	my $rp_c = 0;
	PARA: foreach my $cell ( @chars ) {
			
		if( $cell eq "(" ) {
			push(@bracket_order, "+");
			$lp_c++;
		}
			
		if( $cell eq ")" ) {
			push(@bracket_order, "-");
			$rp_c++;
		}
			
		last PARA if $lp_c == $rp_c;
	}
	
	$n_count = $lb_c + 1;
	$l_count = ( 2 * $lb_c ) + 1;
	
	# collect numerical values
	my @entries = split(/[\(\)\:\,;]+/, $line);
	
	# remove extraneous starting base
	shift(@entries);

	# separate into root and main entries
	my $n_c = 0;
	my $l_c = 0;
	for( my $i = 0; $i < scalar(@entries); $i++) {
		
		if ( $entries[$i] =~ /\./ && $l_c < $l_count ) {
			push(@root_l, $entries[$i]);
			$l_c++;
		} elsif ( $entries[$i] =~ /\./ && $l_c >= $l_count ) {
			push(@main_l, $entries[$i]);
			$l_c++;
		} elsif ( $entries[$i] !~ /\./ && $n_c < $n_count ) {
			push(@root_n, $entries[$i]);
			$n_c++;
		} elsif ( $entries[$i] !~ /\./ && $n_c >= $n_count ) {
			push(@main_n, $entries[$i]);
			$n_c++;
		}
		
	}

	# make long-term variables
	@root_n_persistent = @root_n;
	@main_n_persistent = @main_n;
	
	# count branch lengths
	my $total_root = 0;
	my $total_main = 0;
	( $total_root += $_ ) for @root_l;
	( $total_main += $_ ) for @main_l;
	my $total_length = $total_root + $total_main;
	
	# make lengths array
	my @lengths;
	push(@lengths, @root_l, @main_l);
	
	# calculate TMRCA
	my $tmrca = $lengths[0];
	
	if( @bracket_order ) {
		for( my $x = 1; $x < scalar(@bracket_order); $x++ ) {
		
			if( $bracket_order[($x + 1)] ) {
				if( $bracket_order[$x] eq '-' && $bracket_order[($x+1)] eq '+' ) {
					$tmrca += $lengths[$x];
				}
			}
		}
	$tmrca += $lengths[scalar(@bracket_order)];
	}
	
	# summary variables
	my $root_g1 = 0;
	my $root_g2 = 0;
	my $main_g1 = 0;
	my $main_g2 = 0;
	
	# count root lineages in group one and two
	foreach my $first_one (@root_n) {
		if ( $first_one <= $group_one ) {
			$root_g1++;
		}
		if ( $first_one > $group_one && $first_one <= ($group_one + $group_two) ) {
			$root_g2++;
		}		
	}

	# count main lineages in group one and two
	foreach my $second_one (@main_n) {
		if ( $second_one <= $group_one ) {
			$main_g1++;
		}
		if ( $second_one > $group_one && $second_one <= ($group_one + $group_two) ) {
			$main_g2++;
		}		
	}

	# count root and main lineages
	my $root_c = scalar(@root_n);
	my $main_c = scalar(@main_n);
	
	# calculate proportions
	my $root_prop = $r_prop = sprintf("%.5f", ($root_g1/$root_c));
	my $main_prop = $m_prop = sprintf("%.5f", ($main_g1/$main_c));
	my $max_prop = sprintf("%.5f", &max(($root_g1/$root_c),($main_g1/$main_c)));
	my $min_prop = sprintf("%.5f", &min(($root_g1/$root_c),($main_g1/$main_c)));
	
	# print out summaries
	print STDOUT "$min_prop\t$tmrca\n";
	
}

# min subroutine
sub min($$) {
    if ($_[0]>$_[1]) {return $_[1]} else {return $_[0]};
}

# max subroutine
sub max($$) {
    if ($_[0]<$_[1]) {return $_[1]} else {return $_[0]};
}

