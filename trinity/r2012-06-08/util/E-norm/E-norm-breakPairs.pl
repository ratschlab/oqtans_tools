#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fastq_reader;
use CMD_processor;
use Data::Dumper;

my $usage = "\n\nusage: $0 left.fq right.fq [max_abundance]\n\n";

my $left_fq_file = $ARGV[0] or die $usage;
my $right_fq_file = $ARGV[1] or die $usage;
my $max_abundance = $ARGV[2];

main: {

    &write_read_sorted_tab_files($left_fq_file, $right_fq_file);
 
    my $left_sorted_tab_file = "$left_fq_file.tab.sort";
    my $right_sorted_tab_file = "$right_fq_file.tab.sort";
       
    unless ($max_abundance) {
    
        ## get read count distribution
        print STDERR "-counting read abundance values\n";
        my %read_count_counter;
        &count_read_abundance($left_sorted_tab_file, \%read_count_counter);
        &count_read_abundance($right_sorted_tab_file, \%read_count_counter);
        
        # print Dumper(\%read_count_counter);
        print STDERR "-plotting histogram of counts\n";
        &generate_histogram(\%read_count_counter);
                
        print STDERR "-computing median coverage value\n";
        $max_abundance = &compute_median_val_ignore_singletons(\%read_count_counter);
    
        print STDERR "-setting max abundance to median replication value: $max_abundance\n";
    }
    
    ## write normalized merged tab files.
    print STDERR "-writing normalized tab files\n";
    my ($normalized_left_tab_file, $normalized_right_tab_file) = &write_normalized_files($max_abundance, $left_sorted_tab_file, $right_sorted_tab_file);
    
    print STDERR "-generating final normalized fastq files\n";
    &convert_to_paired_normalized_fq_files($max_abundance, $normalized_left_tab_file, "$left_fq_file.m${max_abundance}.fq");
    &convert_to_paired_normalized_fq_files($max_abundance, $normalized_right_tab_file, "$right_fq_file.m${max_abundance}.fq");
        
    print STDERR "-done.\n\n";
    
    exit(0);

}


####
sub generate_histogram {
    my ($read_count_counter_href) = @_;

    open (my $ofh, ">counts.hist") or die $!;
    foreach my $bin (sort {$a<=>$b} keys %$read_count_counter_href) {
        
        my $val = $read_count_counter_href->{$bin};
        print $ofh "$bin\t$val\n";
    }
    close $ofh;

    open ($ofh, ">counts.hist.R") or die $!;
    print $ofh "data = read.table(\"counts.hist\", header=F)\n";
    print $ofh "pdf(\"counts.hist.pdf\")\n";
    print $ofh "plot(data[,1], log2(data[,2]), xlab='read abundance bin', ylab='log2(count of reads)')\n";
    print $ofh "dev.off()\n";
    close $ofh;
    
    my $cmd = "R --vanilla -q < counts.hist.R";
    &process_cmd($cmd);
    
    return;
}




####
sub convert_to_paired_normalized_fq_files {
    my ($cutoff_val, $normalized_tab_file, $output_fq_filename) = @_;

    open (my $fh, $normalized_tab_file) or die "Error, cannot open file $normalized_tab_file";

    open (my $ofh, ">$output_fq_filename") or die "Error, cannot write to $output_fq_filename";
    
    while (<$fh>) {
        chomp;
        my ($header, $seq, $qual) = split(/\t/);
        $header =~ s/$;/ /g;
        
        if ($seq ne ".") {
            print $ofh join ("\n", $header, $seq, "+", $qual) . "\n";
        }
    }
    close $fh;
    close $ofh;
    
    return;
    
}



####
sub count_read_abundance {
    my ($sorted_tab_file, $read_count_counter_href) = @_;

    my $read_index = 1;

    my $prev_seq = "";
    my $prev_count = 0;
    
    open (my $fh, $sorted_tab_file) or die "Error, cannot open file $sorted_tab_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $read_seq = $x[$read_index];
        
        if ($read_seq eq ".") { next; }
        
        if ($prev_seq eq $read_seq) {
            $prev_count++;
        }
        else {
            if ($prev_count > 0) {
                $read_count_counter_href->{$prev_count}++;
            }
            $prev_seq = $read_seq;
            $prev_count = 1;
        }
    }
    
    # get last one
    $read_count_counter_href->{$prev_count}++;
    
    close $fh;
 
    return;
}


####
sub compute_median_val_ignore_singletons {
    my ($read_count_counter_href) = @_;
    
    my @bins = sort {$a<=>$b} keys %$read_count_counter_href;    
  
    my $sum_counts = 0;
    foreach my $bin (@bins) {
        if ($bin == 1) { next; }
        my $val = $read_count_counter_href->{$bin};
        $sum_counts += $bin * $val;
    }
    
    my $mid_pos = $sum_counts/2;

    print "mid pos = $mid_pos\n";

    my $running_total = 0;
    foreach my $bin (@bins) {

        if ($bin == 1) { next; }

        my $val = $read_count_counter_href->{$bin};
        $running_total += $bin * $val;
        
        print "$bin\t$running_total\t looking for $mid_pos\n";
        
        
        if ($mid_pos <= $running_total) {
            return($bin);
        }
        
        
    }
    
    print STDERR "** warning **, no min count computed: setting to 1.\n";
    
    return(1);
}


####
sub sort_merged_tab_by_read_seqs {
    my ($merged_tab_file, $left_read_sorted_merged_file, $right_read_sorted_merged_file) = @_;

    my @cmds;
    push (@cmds, "sort -T . -S 2G -k 2,2 $merged_tab_file > $left_read_sorted_merged_file") unless (-s "$left_read_sorted_merged_file" && -s "$left_read_sorted_merged_file" == -s "$merged_tab_file");
    
    push (@cmds, "sort -T . -S 2G -k 5,5 $merged_tab_file > $right_read_sorted_merged_file") unless (-s "$right_read_sorted_merged_file" && -s "$right_read_sorted_merged_file" == -s "$merged_tab_file");
    
    &process_parallel_cmds(@cmds) if @cmds;
    
    return;
}


####
sub write_read_sorted_tab_files {
    my ($left_fq_file, $right_fq_file) = @_;
    

    ########################################################
    ## Generate individual sorted tab files for each fq file
    ########################################################


    ## ok, here's where I could use some multi-threading instead.
    
    foreach my $fq_file ($left_fq_file, $right_fq_file) {
        my $pid = fork();
        unless ($pid) {
            ## child
            &make_read_sorted_tab_file($fq_file);
            
            exit(0);
        }
    }
    
    my $something_failed = 0;
    
    while (my $pid = wait()) {
        if ($pid > 0) {
            if ($?) {
                print STDERR "Error, process $pid failed with ret $?\n";
                $something_failed++;
            }
        }
        else {
            last;
        }
    }
    
    if ($something_failed) {
        die "Process failed somewhere. See error msgs above.\n";
    }
        
    return;
        
}



####
sub make_read_sorted_tab_file {
    my ($fq_file) = @_;

    
    my $fq_tab_filename = "$fq_file.tab";


    unless (-e "$fq_tab_filename.ok") {
    
        print STDERR "-writing tab file: $fq_tab_filename\n";
        
        open (my $ofh, ">$fq_tab_filename") or die "Error, cannot write to $fq_tab_filename";
        
        my $fq_reader = new Fastq_reader($fq_file);
        while (my $record = $fq_reader->next()) {
            my ($acc, $seq, $plus, $qual) = split(/\n/, $record->get_fastq_record());
            $acc =~ s/\s/$;/g;
            print $ofh join("\t", $acc, $seq, $qual) . "\n";
        }
        close $ofh;
        
        &process_cmd("touch $fq_tab_filename.ok");
    }

    my $sorted_tab_file = "$fq_tab_filename.sort";
    unless (-e "$sorted_tab_file.ok") {
        ## sort by read sequence
        my $cmd = "sort -T . -S 2G -k2,2 $fq_tab_filename > $sorted_tab_file";
        &process_cmd($cmd);
        
        &process_cmd("touch $sorted_tab_file.ok");
    }
    
    return;
    
}


####
sub write_normalized_files {
    my ($median_val, $left_read_sorted_tab_file, $right_read_sorted_tab_file) = @_;

    
    my $left_normalized_tab_file = "$left_read_sorted_tab_file.m$median_val.norm";
    my $right_normalized_tab_file = "$right_read_sorted_tab_file.m$median_val.norm";
    

    my $pid = fork();
    unless ($pid) {
        &write_normalized_file($median_val, $left_read_sorted_tab_file, $left_normalized_tab_file) unless (-s $left_normalized_tab_file);
        exit(0);
    }

    $pid = fork();
    unless ($pid) {
        &write_normalized_file($median_val, $right_read_sorted_tab_file, $right_normalized_tab_file) unless (-s $right_normalized_tab_file);
        exit(0);
    }
    
    # cheap multithreading
    
    my $failed = 0;
    while ($pid = wait()) {
        if ($pid <= 0) {
            last;
        }
        else {
            if ($?) {
                print STDERR "Error, process $pid failed with ret $?\n";
                $failed = 1;
            }
        }
    }
    if ($failed) {
        die "Error, at least one process failed, see above.";
    }
    
    return ($left_normalized_tab_file, $right_normalized_tab_file);
}

####
sub write_normalized_file {
    my ($max_cutoff, $filename, $outputfilename) = @_;
    
    my $read_index = 1; 
    
    my $prev_seq = 0;
    my $prev_count = 0;
    
    open (my $ofh, ">$outputfilename") or die "Error, cannot write to $outputfilename";
    open (my $fh, $filename) or die "Error, cannot open file $filename";
    while (<$fh>) {
        my $line = $_;
        chomp;
        my @x = split(/\t/);
        my $read_seq = $x[$read_index];
        if ($prev_seq eq $read_seq) {
            $prev_count++;
            if ($prev_count <= $max_cutoff) {
                print $ofh $line;
            }
        }
        else {
            $prev_seq = $read_seq;
            $prev_count = 1;
            print $ofh $line;
        }
    }
    close $fh;
    close $ofh;
    
    return;
}



####
package Tab_reader;

use strict;
use warnings;

sub new {
    my ($packagename) = shift;
    my ($filename) = @_;
        
    my $self = {
        filename => $filename,
        fh => undef,
    };
    bless ($self, $packagename);
    
    open (my $fh, $filename) or die $!;
    $self->{fh} = $fh;
    
    return($self);
}

sub get_next {
    my $self = shift;

    my $fh = $self->{fh};
    my $line = <$fh>;
        
    if ($line) {
        chomp $line;
        my $tab_obj = Tab->new($line);
        return($tab_obj);
    }
    else {
        return(undef);
    }
}

############################
package Tab;

use strict;
use warnings;

sub new {
    my $packagename = shift;
    my ($line) = @_;

    my ($header, $seq, $qual) = split(/\t/, $line);
    $header =~ s/$;/ /g;
    
    my $self = { header => $header,
                 seq => $seq,
                 qual => $qual,
             };

    bless ($self, $packagename);
    
    return($self);
}


####
sub get_header {
    my $self = shift;

    return($self->{header});
}

####
sub get_full_accession {
    my $self = shift;

    my $header = $self->get_header();

    my @fields = split(/\s+/, $header);
    return($fields[0]);
}

####
sub get_core_accession {
    my $self = shift;

    my $accession = $self->get_full_accession();
    $accession =~ s/\/[12]$//;
    
    return($accession);

}

####
sub tab_entry {
    my $self = shift;
    
    my $header = $self->get_header();
    $header =~ s/\s/$;/g;
    
    return(join("\t", $header, $self->{seq}, $self->{qual}));
}
