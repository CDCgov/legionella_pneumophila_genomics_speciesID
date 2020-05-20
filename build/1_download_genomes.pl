use strict;
use warnings;
use LWP::Simple;
use Cwd;

#  This script will get all Legionella sp genomes between 500k and 7m in length:

#	Plasmid specific sequences will be avoided, but I make no promise that the resulting
#	set will be plasmid free.

#	Using code found here as a starting point from NCBI:
#	https://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.Application_3_Retrieving_large

#	For ENA, pull,

#  ENA taxon ID for Legionella genus
my $ena_taxid = 445;

my $out_dir = "ref_genomes";


#  To deal with some quirks in how data is returned, I had to adjust
#  my original script.  Rather than change my data structure
#  to deal with multiple contigs, I concatenate to the fasta files I 
#  create.  If there are already files there, it could make a big mess.
#  Thus I push this off on the user.  They should back up their old
#  reference files anyways.  :-)

#  Further note: Merging contigs into the right files doesn't work anyways
#  There is no standard format for naming contigs/scaffolds from same source
#  and the end result would not improve MASH results anyways.

#  UPDATE: This potentially could improve MASH results.  Testing shows that
#  if "contig 1" is searched, it doesn't match "contig 2" in the database
#  as well as it matches complete sequences from some other species.  In a
#  real shotgun sequencing experiment rather than a test, this is an extremely
#  unlikely event, however I will try to clean up the input data for DB 
#  generation to eliminate this as a possibility.

#  UPDATE 2: ENA sequences are all redundant with NCBI seqs.  ENA sequence
#  download is disabled in current version
 


if (-e $out_dir) {

	die ("\n\n\"$out_dir\" directory already exists.\n",
		"This script will not run with that directory present.\n",
		"Please delete or rename \"$out_dir\" and restart script\n\n\n");
} 

mkdir $out_dir;

my $cwd = getcwd();

$out_dir = $cwd . "/" . $out_dir;

#  I'm going to hold all query results seqs in memory before printing so I 
#  can make sure I'm not clobbering names or duplicating.
my $all_seqs;

#  Fallthrough logic on ENA scaffolds is cleaner if I store them as a separate thing.
my $alt_seqs;

my $overwrite_counter = 0;

#  Count species retrieved and sequence keep names
my $species_counts;

#  Keep track of all stripped ID #s to avoid duplication
#  I have observed that the same entry can be downloaded
#  with both an "NZ_" prefix and without.  Track ids
#  so we can skip duplicates.
my %seen_ids;


#  All Legionella sp. sequences with 1,000,000-7,000,000 bases total length
my $query = "Legionella[orgn]+AND+1000000:7000000[SLEN]";

DoQueryNCBI2 ($query);

$query = "Fluoribacter[orgn]+AND+1000000:7000000[SLEN]";

DoQueryNCBI2 ($query);

$query = "Tatlockia[orgn]+AND+1000000:7000000[SLEN]";

DoQueryNCBI2 ($query);


# DoQueryENA ($ena_taxid);

PrintFastaFiles($all_seqs, $out_dir, 1);

# PrintAltFastaFiles($alt_seqs, $out_dir);


my $dt = GetDateTime();
my $stat_file = "reference_retrieval_stats_" . $dt -> [0] . ".txt";
open (my $statfh, ">", $stat_file) or die $!;
foreach my $species (sort keys %$species_counts) {
	my $count = @{$species_counts -> {$species}};
	print $statfh "$species\t$count\n";
	foreach my $strain (sort @{$species_counts -> {$species}}) {
		print $statfh "\t$strain\n";
	}
}
	

sub DoQueryENA {
	my $taxid = shift;
	
	#  Rather than parse their dreadful XML, we're just going to do it "wrong"
	#  and request "text" results.  Then we'll parse the errors reported
	my $url = "https://www.ebi.ac.uk/ena/data/view/Taxon:$taxid&result=assembly&subtree=true&display=text";
	
	my $result = get($url);
	
	my @rows = split ("\n", $result);
	my @ids;
	foreach my $row (@rows) {
		$row =~ s/\R+//g;
		$row =~ m/Display txt not supported for entry (\S+)/;
		push (@ids, $1);
	}
	my $counter = 0;
	foreach my $id (@ids) {
		print "Retrieving $id\n";
		my $url = "https://www.ebi.ac.uk/ena/data/view/$id&display=xml";
		my $result = get ($url);
# 		print "$result\n";
		my @chroms = ($result =~ m/<CHROMOSOMES>\s+<CHROMOSOME accession=\"(\S+)\">\n\s*<TYPE>Chromosome<\/TYPE>/gi);
	
		#  0 chromosomes means fall through to ENA sequence download
		if (@chroms == 0) {
			EnaAltDownload ($result);
			$counter++;		
			next;
		}
		if (@chroms > 1 ) {
			print "\tToo many chromosomes identified\n";
			die;
		}		
		foreach my $chrom (@chroms) {
# 		my $chrom = $1;
			if (defined $chrom) {
				print "\t$chrom\n";
				
				my $url = "http://www.ebi.ac.uk/ena/data/view/$chrom\&display=fasta";
				my $result = get ($url);
							
		        my $seqs = ParseFastaString ($result);
		        foreach my $seq (keys %$seqs) {
			        my $new_name = $seq;
			        $new_name =~ m/(.*)\|(.*)\|(\S+)/;
			        my $id = $2;
			        $new_name =~ s/(.*)\|(.*)\|(\S+)/$1_$3/gi;
					print "\t$id\t$seq\t$new_name\n";

					foreach my $key (keys %$all_seqs) {
						if ($key =~ m/$id/gi) {
							print "Identified $new_name as an ENA duplicate of NCBI sequence $key\n";
							next;
						}
					}
			        if (exists $all_seqs -> {$new_name}) {
				        print "Overwriting $new_name in memory\n";
				        $overwrite_counter++;
			        }
			        $all_seqs -> {$new_name} = $seqs -> {$seq};  				
				
				}
			}
			else {
				print "$result\n";
			}
		}
		$counter++;
		if ($counter > 20) {
			last;
		}
		
	}
# 	my $count = @ids;
# 	print "$count\n";

}



sub ParseFastaString {
    my $string = shift;
    my @lines = split ("\n", $string);
    #print "$file\n";
    my $align;
    my $otu = "";
    foreach my $line (@lines) {
        chomp $line;
        $line =~ s/\R//g;
        if ($line =~ m/^>(.+)/) {
            $otu = $1;
        }
        else {
            #This weird thing is because MACSE indicates frameshifts with a "!".  
            #EVERYTHING chokes on a "!" in a sequence!
            #Replace "!" with gaps
            $line =~ s/\!/-/gi;
            $align -> {$otu} .= $line;
        }
    }
    return $align;
} 

sub PrintFastaFiles {
	my $seqs = shift;
	my $out_dir = shift;
	my $attempt = shift;
	my $filename = shift;	#Should only be passed on attempt=2.  I realize code is getting pretty spaghettified....

	foreach my $key (sort keys %$seqs) {
# 		print "\n!$key\t$attempt\n";
# 		$key =~ s/Tatlockia/Legionella/gi;
# 		$key =~ s/Fluorobacter/Legionella/gi;		
		unless (($key =~ m/^(\S+\.\d+) (Legionella \S+)/i) || ($key =~ m/^(\S+\.\d+) Tatlockia/gi) || ($key =~ m/^(\S+\.\d+) Fluoribacter/gi)) {
			if ($attempt <= 1) {
				print "\t\"$key\" does not match name parsing regex.  Skipping.\n";
				next;
			}
		}
		
		if ($key =~ m/plasmid/i) {
			print "\t\"$key\" appears to be a plasmid.  Skipping.\n";
			next;
		}
		
		if (($key =~ m/Legionella sp.{0,1}$/i) || ($key =~ m/Legionella genomosp.{0,1}$/i)) {
			print "\t\"$key\" corresponds to an unknown Legionella species.  Skipping.\n";
			next;
		}
		
		#Process $key and make filename
		my $acc;
		my $species;		
		if ($key =~ m/^(\S+\.\d+)\s*(Legionella \S+)/i){
			$acc = $1;
			$species = $2;
		}
		elsif ($key =~ m/^(\S+\.\d+)\s*(Fluoribacter \S+)/i){
			$acc = $1;
			$species = $2;
# 			$species =~ s/Fluoribacter/Legionella/gi;
		}
		elsif ($key =~ m/^(\S+\.\d+)\s*(Tatlockia \S+)/i){
			$acc = $1;
			$species = $2;
# 			$species =~ s/Tatlockia/Legionella/gi;
		}
		elsif ($key =~ m/^(\S+\.\d+)\s*(.*)/i){
			$acc = $1;
			$species = $2;
		}
		
		$species = FixBinomialName($species);
			
		unless ( (defined $acc) && (defined $species)) {
			print "regex fail for $key\n";
		}			
		
# 		my $strain = $3;
# 		my $type = $4;

		#  Print friendly versions of fields		
		my $pacc = $acc;
		$pacc =~ s/ /_/g;
		my $pspecies = $species;
		$pspecies =~ s/ /_/g;
		$pspecies =~ s/,//g;
		$pspecies =~ s/\.//g;		
# 		my $pstrain = $strain;
# 		$strain =~ s/ /_/g;
# 		my $ptype = $type;
# 		$ptype =~ s/ /_/g;

		unless ($attempt > 1) {
			$filename = $out_dir . "/" . $pspecies . "-$pacc" . ".fas";
			$pacc =~ m/([^_]*)$/;
			my $id = $1;
			if (exists $seen_ids{$id}) {
				next;
			}
			else {
				$seen_ids{$id} = 1;
			}
		}						
		
		unless (defined $seqs -> {$key}) {
			print "\"$key\" is undefined. Skipping sequence\n";
			next;
		}
		
		#Skip empty contigs
		unless ($seqs -> {$key} =~ m/[^N]/i) {
			if ($attempt > 1) {
				print "\"$key\" contains no sequence data. Skipping sequence\n";
				next;
			}
				
			print "\"$key\" contains no sequence data.  Pushing to secondary download strategy\n";
			my $seqs = GetContigsFromMaster($key);
			PrintFastaFiles($seqs, $out_dir, "2", $filename);
			#  Track what we've seen
			push (@{$species_counts -> {$pspecies}}, $key);	
			next;			
		}			
		
	
		#  Track what we've seen
		if ($attempt <2) {
			push (@{$species_counts -> {$pspecies}}, $key);
		}
					
# 		print "\t $filename\n";
		open (my $out_fh, ">>", $filename) or die $!;
		print $out_fh ">$key\n", $seqs -> {$key}, "\n"; 
	}
}
	
sub GetDateTime {
	#  The Perl module Date::Time makes much more sense to use, but it isn't
	#  installed on all of our systems.  This will work on anything (until it doesn't...)
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime();

	$mon += 1; #Month returned in 0-11 range
	$year += 1900; #Year returned in years since 1900
	
	#Pad with zeros
	$mon = sprintf "%02s", $mon;
	$mday = sprintf "%02s", $mday;
	$hour = sprintf "%02s", $hour;
	$min = sprintf "%02s", $min;
	$sec = sprintf "%02s", $sec;
	
	my $date = "$year-$mon-$mday";
	my $time = "$hour:$min:$sec";
	my @dt = ($date,$time);
	
	return \@dt;
}
sub DoQueryNCBI2 {
	my $query = shift;
	
	my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
	my $url = $base . "esearch.fcgi?db=nucleotide&term=$query&retmax=100000&usehistory=y";
	
	my $output = get($url);
	
	#parse WebEnv, QueryKey and Count (# records retrieved)
	my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
	my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
	my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);
	
	# print "$output\n\n";
	
# 	print "$web\t$key\t$count\n";
	
	#retrieve data in batches of 500
	
	
	my $retmax = 500;
	for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
        my $efetch_url = $base ."efetch.fcgi?db=nucleotide&WebEnv=$web";
        $efetch_url .= "&query_key=$key&retstart=$retstart";
        $efetch_url .= "&retmax=$retmax&rettype=fasta&retmode=text";
        my $efetch_out = get($efetch_url);
#         print "$efetch_out\n";

        my $seqs = ParseFastaString ($efetch_out);
        	
        foreach my $seq (keys %$seqs) {
	        if ($seq =~ m/[^\d]0[1-9]000\d+/) {
	        	print "Skipping $seq - Identified as likely contig associated with another file\n";		        
		        next;
	        }
	      	unless (($seq =~ m/complete genome/) || ($seq =~ m/chromosome/) || ($seq =~ m/whole genome/)) {
	        	print "Skipping $seq\n";		        
		        next;
	        }
	        #Skip unassigned species
	        if ($seq =~ m/Legionella sp\. /) {
	        	print "Skipping $seq\n";		        
		        next;
	        }
	        #Same unassigned Legionella sauodensis as above.  Tough bug to skip!
	        if ($seq =~ m/LN901324\.1/) {
	        	print "Skipping $seq\n";		        
		        next;
	        } 	         
	        
	        
# 	        unless ($seq =~ m/00000000/i) {
# 		        unless ($seq =~ m/complete genome/i) {
# 		        	print "Skipping $seq\n";
# 		        	next;
# 	        	}
# 	        }
		        	
	        if (exists $all_seqs -> {$seq}) {
		        print "Overwriting $seq in memory\n";
		        $overwrite_counter++;
	        }
	        $all_seqs -> {$seq} = $seqs -> {$seq};     
		}
# 		last;
	}
}

sub GetContigsFromMaster {
	my $m_id = shift;
	print "$m_id\n";
	
	#  Searching on whole sequence name is not always working for inscrutable
	#  reasons.  Extract just the ID.
	$m_id = FixBinomialName($m_id);
	my $acc;
	if ($m_id =~ m/^(\S+)\s+Legio/i) {
		$acc = $1;
	}
	elsif ($m_id =~ m/^(\S+)\s+Fluori/i) {
		$acc = $1;
	}
	elsif ($m_id =~ m/^(\S+)\s+Tatlock/i) {
		$acc = $1;
	}		
	my $query = $acc;
	my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
	my $url = $base . "esearch.fcgi?db=nucleotide&term=$query&retmax=100000&usehistory=y";

	my $try = 0;
	my $web;
	my $key;
	my $count;
	while ($try < 5) {
		my $output = get($url);	
			
		$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
		$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
		$count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);
		
		if ((defined $web) && (defined $key) && (defined $count)) {
			last;
		}
		print STDERR "$output\n";
		sleep 5;
		$try++
	}	
	if ($try == 5) {
		die ("Sequence for $m_id could not be determined");
	}
	
    if ($count > 1) {
	    die ("More than 1 result returned for $m_id\n");
    }
#     print "$web\t$key\t$count\n";   

    #Get GB file for master WGS page	
    my $efetch_url = $base ."efetch.fcgi?db=nucleotide&WebEnv=$web";
    $efetch_url .= "&query_key=$key&rettype=gb&retmode=text";
    my $efetch_out = get($efetch_url);
#     print "$efetch_out\n";
    
    #Parse out range of associated contigs
    unless ($efetch_out =~m/^WGS\s+(\S+)-(\S+)/m){
	    print "$m_id gb does not contain contig ID range\n";
	    return;
	    
    }     
    $efetch_out =~m/^WGS\s+(\S+)-(\S+)/;
    
    my $start = $1;
    my $end = $2;
    
    #Construct new query
    $query = "$start:$end" . "[PACC]";
    
	$url = $base . "esearch.fcgi?db=nucleotide&term=$query&retmax=100000&usehistory=y";
	$try = 0;
	$web = '';
	$key = '';
	$count = '';
	while ($try < 5) {
		my $output = get($url);	
			
		$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
		$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
		$count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);
		
		if ((defined $web) && (defined $key) && (defined $count)) {
			last;
		}
		print STDERR "$output\n";
		sleep 5;
		$try++
	}	
	if ($try == 5) {
		die ("Sequence for $m_id could not be determined");
	}
	
	my $retmax = 500;
	for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
        my $efetch_url = $base ."efetch.fcgi?db=nucleotide&WebEnv=$web";
        $efetch_url .= "&query_key=$key&retstart=$retstart";
        $efetch_url .= "&retmax=$retmax&rettype=fasta&retmode=text";
        my $efetch_out = get($efetch_url);

        $efetch_out = FixBinomialName($efetch_out);
        my $seqs = ParseFastaString ($efetch_out);
		return $seqs;        

	}
}

sub EnaAltDownload {
	my $result = shift;
	print "\tNo chromosomes identified\n";
	$result =~ m/<TAXON_ID>(\d+)<\/TAXON_ID>/gi;
	my $iso_id = $1;
	print "\tBegining alternate download strategy for ENA sequences associated with TAXON ID $iso_id\n";

	my $primary_id = "Taxid_$iso_id";
	if ($result =~ m/<PRIMARY_ID>(\S+)<\/PRIMARY_ID>/gi) {
		$primary_id = $1;
	}
	elsif ($result =~ m/accession=\"(\S+)\"/gi) {
		$primary_id = $1;
	}
	my $title = "XXX";
	if ($result =~ m/<TITLE>([^\n]+)<\/TITLE>/gi) {
		$title = $1;
	}
	elsif ($result =~ m/<SCIENTIFIC_NAME>([^\n]+)<\/SCIENTIFIC_NAME>/gi) {
		$title = $1;
	}
	
	my $new_name = $primary_id . " " . $title;
	print "\tNew name: $new_name\n";	
	
	my $url = "http://www.ebi.ac.uk/ena/data/view/Taxon:$iso_id&result=sequence_release&display=fasta";
	my $new_result = get ($url);
	
	#  No result
	if ($new_result =~ m/display type is either not supported or entry is not found/i) {
		print "\tEntry not found\n";
		return;
	}
	
	#  Crude check for adequate length
	if (length ($new_result) < 1000000) {
		print "\tSequence too short\n";
		return;
	}
		
	
	my $seqs = ParseFastaString ($new_result);
	
	$alt_seqs -> {$new_name} = $seqs;
}
sub PrintAltFastaFiles {
	my $seqs = shift;
	my $out_dir = shift;

	foreach my $key (sort keys %$seqs) {
# 		print "\n!$key\t$attempt\n";
		unless ($key =~ m/(\S+).*(Legionella \S+)/i) {
			print "\t\"$key\" does not match name parsing regex.  Skipping.\n";
			next;
		}
		
		if ($key =~ m/plasmid/i) {
			print "\t\"$key\" appears to be a plasmid.  Skipping.\n";
			next;
		}
		
		#Process $key and make filename
		$key =~ m/(\S+).*(Legionella \S+)/i;
		my $acc = $1;
		my $species = $2;
# 		my $strain = $3;
# 		my $type = $4;

		#  Print friendly versions of fields		
		my $pacc = $acc;
		$pacc =~ s/ /_/g;
		my $pspecies = $species;
		$pspecies =~ s/ /_/g;
		$pspecies =~ s/,//g;
		$pspecies =~ s/\.//g;		
# 		my $pstrain = $strain;
# 		$strain =~ s/ /_/g;
# 		my $ptype = $type;
# 		$ptype =~ s/ /_/g;

		my $filename = $out_dir . "/" . $pspecies . "-$pacc" . ".fas";						
		open (my $out_fh, ">>", $filename) or die $!;
		foreach my $contig (keys %{$seqs -> {$key}}) {
			print $out_fh ">$contig\n", $seqs -> {$key} -> {$contig}, "\n";
		} 		
		#  Track what we've seen

		push (@{$species_counts -> {$pspecies}}, $key);

					
# 		print "\t $filename\n";
	}		
}	

sub FixBinomialName {
	my $species = shift;
	$species =~ s/Fluoribacter/Legionella/gi;
	$species =~ s/Tatlockia/Legionella/gi;
	unless ($species =~ m/bozemanii/i) {
		$species =~ s/bozemani/bozemanii/gi;	
		$species =~ s/bozemanae/bozemanii/gi;
	}
	$species =~ s/Legionella subs pascullei/Legionella pneumophila/gi;			
	$species =~ s/Legionella_subs_pascullei/Legionella_pneumophila/gi;
	
	return $species;
}	
