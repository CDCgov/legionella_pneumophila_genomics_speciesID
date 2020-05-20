use strict;
use warnings;
use File::Find;
use Cwd;

#  Run this script on cluster after downloading all genomes to /ref_genomes
#  Generates individual MASH sketches for each genome in /mash_sketches

my $kmer_size = 25;
my $sketch_size = 100000;
my $threads = 8;

my @source_dirs = qw (ref_genomes ref_genomes_rare);
my $out_dir = "mash_sketches";

my $cwd = getcwd();
$out_dir = $cwd . "/" . $out_dir;

if (-e $out_dir) {
	die ("\n\n\"$out_dir\" directory already exists.\n",
		"This script will not run with that directory present.\n",
		"Please delete or rename \"$out_dir\" and restart script\n\n\n");
} 

# Get fasta files
my @fasta_files;
foreach my $source_dir (@source_dirs) {
	my $dir = $cwd . "/" . $source_dir;
	chdir $dir;
	print "$dir\n";
	my @files = <*.fas>;
	@files = map {$dir . "/" . $_} @files;
	push (@fasta_files, @files);
}
mkdir $out_dir;	 
chdir $out_dir;

foreach my $file (@fasta_files) {
	
	$file =~ m/([^\/]*)$/;
	my $file_root = $1;
	print "$file_root\n";
	my $full_fasta_path = $file;
	my $full_output_path = $out_dir . "/" . $file_root;
	
# 	my $command_string = "mash sketch -k $kmer_size -s $sketch_size -p $threads";
# 	$command_string .= " -o $full_output_path $full_fasta_path";
# 	my $result = `$command_string`;
	
	$file_root =~ s/\.fas//;
	
	MakeSGESub ($file_root, $full_fasta_path, $full_output_path);
# 	sleep (1);
# 	last;

	
}

sub MakeSGESub {
    my $file_root = shift;
    my $full_fasta_path = shift;
    my $full_output_path = shift;

    my $type = "mash";
    my $threads = 8;
    
    my $script_file = $file_root . "_$type.sge";
    
    open (my $out_fh, ">", $script_file) or die $!;
    print $out_fh "\#!\/bin\/bash -l\n";
    print $out_fh "\n";
    print $out_fh "\# Embedded Grid Engine parameters\n";
    print $out_fh "\# All Embedded parameters must start with \"#$\"\n";
    print $out_fh "\# This is the same as adding these lines to the actual qsub line\n";
    print $out_fh "\n";
    print $out_fh "\# save the standard output text to this file instead of the  the default jobID.o file\n";
    print $out_fh "\#\$ -o $file_root", "_$type.out\n";
    print $out_fh "\n";
    print $out_fh "\# save the standard error text to this file instead of the  the default jobID.e file\n";
    print $out_fh "\#\$ -e $file_root", "_$type.err\n";
    print $out_fh "\n";
    print $out_fh "\# Rename the job to be this string instead of the default which is the name of the script\n";    
    print $out_fh "\#\$ -N $file_root", "_$type\n";
    print $out_fh "\n";
    print $out_fh "\# Threaded job $threads slots on the same node\n";
    print $out_fh "\#\$ -pe smp $threads\n";
    print $out_fh "\n";        
    print $out_fh "\# Refer all file reference to work the current working directory which is\n";
    print $out_fh "\# the directory from which the script was qsubbed\n";
    print $out_fh "\#\$ -cwd\n";
    print $out_fh "\n";
    print $out_fh "\# Always run in the default queue\n";
    print $out_fh "\#\$ -q all.q\n";       
    print $out_fh "\n";
    print $out_fh "\# Load the modules\n";
    print $out_fh "module load Mash/2.0\n";
            
    print $out_fh "\n";
    print $out_fh "\# perform the actual $type command with appropriate command line options\n";
    print $out_fh "\n";
    
	my $command_string = "-k $kmer_size -s $sketch_size -p $threads";
	$command_string .= " -o $full_output_path $full_fasta_path";    

	print $out_fh "mash sketch $command_string\n";
	print $out_fh "\n";

    print $out_fh "\# unload the module now that we are done using the commands\n";
    print $out_fh "module unload Mash/2.0\n";             
    
    close $out_fh;
    `qsub $script_file`;
}
 
	







