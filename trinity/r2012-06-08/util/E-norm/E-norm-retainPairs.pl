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

    my $merged_tab_file = "fqtab.merged";
    print STDERR "-writing merged tab file\n";
    &write_merged_tab_file($left_fq_file, $right_fq_file, $merged_tab_file) unless (-s $merged_tab_file);

    ## sort each file by sequence
    print STDERR "-sorting each tab file by sequence\n";
    my $left_read_sorted_merged_file = "$merged_tab_file.left.sort";
    my $right_read_sorted_merged_file = "$merged_tab_file.right.sort";
    &sort_merged_tab_by_read_seqs($merged_tab_file, $left_read_sorted_merged_file, $right_read_sorted_merged_file) unless (-s $left_read_sorted_merged_file && -s $right_read_sorted_merged_file);
    
    
    unless ($max_abundance) {
    
        ## get read count distribution
        print STDERR "-counting read abundance values\n";
        my %read_count_counter;
        &count_read_abundance($left_read_sorted_merged_file, "left", \%read_count_counter);
        &count_read_abundance($right_read_sorted_merged_file, "right", \%read_count_counter);
        
        # print Dumper(\%read_count_counter);
        print STDERR "-plotting histogram of counts\n";
        &generate_histogram(\%read_count_counter);
        
        
        print STDERR "-computing median coverage value\n";
        $max_abundance = &compute_median_val_ignore_singletons(\%read_count_counter);
    
        print STDERR "-setting max abundance to median replication value: $max_abundance\n";
    }
    
    ## write normalized merged tab files.
    print STDERR "-writing normalized tab files\n";
    my $normalized_left_merged_file = "$left_read_sorted_merged_file.norm";
    my $normalized_right_merged_file = "$right_read_sorted_merged_file.norm";
    
    &write_normalized_files($max_abundance, $left_read_sorted_merged_file, $right_read_sorted_merged_file) unless (-s $normalized_left_merged_file && -s $normalized_right_merged_file);
    
    ## merge and unique the results.
    print STDERR "-merging normalized files\n";
    my $cmd = "cat $left_read_sorted_merged_file $right_read_sorted_merged_file | sort -T . -S 2G -k1,1 | uniq > $merged_tab_file.norm";
    &process_cmd($cmd);
    
    print STDERR "-generating final normalized fastq files\n";
    &convert_to_paired_normalized_fq_files($max_abundance, $left_fq_file, $right_fq_file, "$merged_tab_file.norm");
    
    
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
    my ($cutoff_val, $left_fq_file, $right_fq_file, $norm_merged_tab_file) = @_;

    open (my $fh, $norm_merged_tab_file) or die "Error, cannot open file $norm_merged_tab_file";
    my $left_norm_fq_file = "$left_fq_file.m$cutoff_val.fq";
    my $right_norm_fq_file = "$right_fq_file.m$cutoff_val.fq";
    
    open (my $left_ofh, ">$left_norm_fq_file") or die "Error, cannot write to $left_norm_fq_file";
    open (my $right_ofh, ">$right_norm_fq_file") or die "Error, cannot write to $right_norm_fq_file";

    while (<$fh>) {
        chomp;
        my ($headerA, $seqA, $qualA, $headerB, $seqB, $qualB) = split(/\t/);
        $headerA =~ s/$;/ /g;
        $headerB =~ s/$;/ /g;

        if ($seqA ne ".") {
            print $left_ofh join ("\n", $headerA, $seqA, "+", $qualA) . "\n";
        }
        if ($seqB ne ".") {
            print $right_ofh join("\n", $headerB, $seqB, "+", $qualB) . "\n";
        }
    }
    close $fh;
    close $left_ofh;
    close $right_ofh;
    
    return;
    
}



####
sub count_read_abundance {
    my ($sorted_merged_file, $read_end_type, $read_count_counter_href) = @_;

    my $read_index = ($read_end_type eq 'left') ? 1 : 4;

    my $prev_seq = "";
    my $prev_count = 0;
    
    open (my $fh, $sorted_merged_file) or die "Error, cannot open file $sorted_merged_file";
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
    
    die "Error, shouldn't get here";
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
sub write_merged_tab_file {
    my ($left_fq_file, $right_fq_file, $merged_tab_file) = @_;
    

    ########################################################
    ## Generate individual sorted tab files for each fq file
    ########################################################


    ## ok, here's where I could use some multi-threading instead.
    
    foreach my $fq_file ($left_fq_file, $right_fq_file) {
        my $pid = fork();
        unless ($pid) {
            ## child
            &make_sorted_tab_file($fq_file);
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


    open (my $ofh, ">$merged_tab_file") or die "Error, cannot write to $merged_tab_file";
    

    my $left_sort_tab_file = "$left_fq_file.tab.sort";
    my $right_sort_tab_file = "$right_fq_file.tab.sort";
            
    my $left_tab_reader = Tab_reader->new($left_sort_tab_file);
    my $right_tab_reader = Tab_reader->new($right_sort_tab_file);

    my $left_tab = $left_tab_reader->get_next();
    my $right_tab = $right_tab_reader->get_next();
    

    while ($left_tab && $right_tab) {

        if ($left_tab->get_core_accession() eq $right_tab->get_core_accession()) {
            
            ## print merged tab entry
            print $ofh join("\t", $left_tab->tab_entry(), $right_tab->tab_entry()) . "\n";
            
            $left_tab = $left_tab_reader->get_next();
            $right_tab = $right_tab_reader->get_next();
        }
        elsif ($left_tab->get_core_accession() lt $right_tab->get_core_accession()) {
            
            print $ofh join("\t", $left_tab->tab_entry(), ".\t.\t.") . "\n";
            $left_tab = $left_tab_reader->get_next();
        }
        elsif ($left_tab->get_core_accession() gt $right_tab->get_core_accession()) {
            print $ofh join("\t", ".\t.\t.", $right_tab->tab_entry()) . "\n";
            $right_tab = $right_tab_reader->get_next();
        }
    }
    
    ## report the remainders
    while ($left_tab = $left_tab_reader->get_next()) {
        print $ofh join("\t", $left_tab->tab_entry(), ".\t.\t.") . "\n";
    }
    while ($right_tab = $right_tab_reader->get_next()) {
        print $ofh join("\t", ".\t.\t.", $right_tab->tab_entry()) . "\n";
    }
    
    close $ofh;
    
    return;
    
}



####
sub make_sorted_tab_file {
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
        my $cmd = "sort -T . -S 2G -k1,1 $fq_tab_filename > $sorted_tab_file";
        &process_cmd($cmd);
        
        &process_cmd("touch $sorted_tab_file.ok");
    }

    return;
    
}


####
sub write_normalized_files {
    my ($median_val, $left_read_sorted_merged_file, $right_read_sorted_merged_file) = @_;


    my $pid = fork();
    unless ($pid) {
        &write_normalized_file($median_val, $left_read_sorted_merged_file, 'left') unless (-s "$left_read_sorted_merged_file.norm");
        exit(0);
    }

    $pid = fork();
    unless ($pid) {
        &write_normalized_file($median_val, $right_read_sorted_merged_file, 'right') unless (-s "$right_read_sorted_merged_file.norm");
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
    
    return;
}

####
sub write_normalized_file {
    my ($max_cutoff, $filename, $read_end_type) = @_;
    
    my $read_index = ($read_end_type eq 'left') ? 1 : 4;
    
    my $prev_seq = 0;
    my $prev_count = 0;
    
    open (my $ofh, ">$filename.norm") or die "Error, cannot write to $filename.norm";
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
