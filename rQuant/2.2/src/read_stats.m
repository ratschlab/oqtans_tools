function read_stats(anno_dir, track, output_file)
% READ_STATS   Generates a statistic about the read alignments and the covered genes.
%
%   read_stats(anno_dir, track, output_file)
%
%   -- input --
%   anno_dir:        directory of genes
%   track:           name of BAM file
%   output_file:     result file with read statistics: number of
%                    reads, median coverage, number of spliced
%                    reads, overlapping annotated introns
% 
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch
%   Copyright (C) 2009-2011 Max Planck Society
%


% rQuant paths
global RQUANT_PATH RQUANT_SRC_PATH

% rQuant version
global RQUANT_VERSION

% SAMTools path
global SAMTOOLS_DIR

addpath(sprintf('%s/mex', RQUANT_PATH));
addpath(sprintf('%s/tools', RQUANT_PATH));
addpath(sprintf('%s', RQUANT_SRC_PATH));

[ret, timedate] = unix('date');
timedate(timedate==sprintf('\n')) = [];
fprintf(1, '\n------------------------------------------------------------- \n');
fprintf(1, ' ReadStats version %s started %s\n', RQUANT_VERSION, timedate) ;
fprintf(1, '------------------------------------------------------------- \n\n');

CFG.samtools_dir = sprintf('%s/', SAMTOOLS_DIR);
CFG.both_strands = 1;
CFG.paired = 0;
CFG.VERBOSE = 1;

%%%%% genome information %%%%%
try
  [s bam_header] = unix(sprintf('%s./samtools view %s -H', CFG.samtools_dir, track));
  assert(s==0);
  fidx1 = strfind(bam_header,'SN'); fidx2 = strfind(bam_header,'LN');
  for c = 1:length(fidx1), 
    contig_names{c} = bam_header(fidx1(c)+3:fidx2(c)-2); 
  end
catch
  fprintf(1, '\ncontig names could not be parsed from BAM file\n');
  return;
end

%%%%% read length, number of mapped reads, bai file %%%%%
CFG.read_len = 0;
for c = 1:length(contig_names),
  fprintf('Checking bam file: contig %i/%i\n', c, length(contig_names));
  CFG.tracks_fn{c} = {track};
  if ~exist(CFG.tracks_fn{c}{1}, 'file'),
    basedir = fileparts(CFG.tracks_fn{c}{1});
    filename = strrep(CFG.tracks_fn{c}{1}, '_files/alignments.bam', '.dat');
    unix(sprintf('mkdir %s; ln -s %s %s',  basedir, filename, CFG.tracks_fn{c}{1}));
  end 
  for f = 1:length(CFG.tracks_fn{c}),
    fname = CFG.tracks_fn{c}{f};
    if ~exist(sprintf('%s.bai', fname), 'file')
      command = sprintf('%s./samtools index %s', CFG.samtools_dir, fname);
      [s m] = unix(command);
      if ~exist(sprintf('%s.bai', fname), 'file')
        if CFG.VERBOSE>0, fprintf(1, '\nbai file could not be created\n'); end
      end
    end
  end
  try
    [read_len mapped_reads(c)] = get_bam_properties(CFG.tracks_fn{c}{1}, CFG.samtools_dir, contig_names{c});
  catch
    read_len = 0;
    mapped_reads(c) = 0;
  end
  CFG.read_len = max(CFG.read_len, read_len);
end


load(sprintf('%s/genes.mat', anno_dir), 'genes');
genes = sanitise_genes(genes, CFG);
for g = 1:length(genes),
  if mod(g, 100)==0, fprintf(1, 'processed %i genes\n', g); end
  genes(g).mean_ec = 0;
  genes(g).introns_covered = false(1,0);
  genes(g).introns_conf = 0;
  try
    [coverage excluded_reads reads_ok introns read_starts] = get_coverage_per_read(CFG, genes(g));  
    % plus strand
    fidx = find(introns(:,4)==0);
    intron_starts{1} = introns(fidx,1); intron_stops{1} = introns(fidx,2); conf{1} = introns(fidx,3);
    % minus strand
    fidx = find(introns(:,4)==1); intron_starts{2} = introns(fidx,1); intron_stops{2} = introns(fidx,2); conf{2} = introns(fidx,3);
  catch
    reads_ok = 0;
  end
  if ~reads_ok,
    fprintf(1, 'coverage could not be loaded for gene %i\n', g);
    continue;
  end
  genes(g).mean_ec = mean(sum(coverage,2));
  intron_list = zeros(0,4); 
  if CFG.both_strands && isfield(genes(g), 'strands') && ~isempty(genes(g).strands)
    strand_str = unique(genes(g).strands);
  else
    strand_str = genes(g).strand;
  end
  assert(length(strand_str)==1|length(strand_str)==2);
  for s = 1:length(strand_str),
    strand = (strand_str(s)=='-') + 1;
    if ~isempty(intron_starts{strand})
      idx = find(intron_starts{strand}>=genes(g).start & intron_stops{strand}<=genes(g).stop);
      if ~isempty(idx)
        intron_list = [intron_list; double([intron_starts{strand}(idx), intron_stops{strand}(idx), ...
                            conf{strand}(idx), strand*ones(length(idx),1)])];
      end
    end
  end
  introns = [];
  for t = 1:length(genes(g).transcripts),
    introns = [introns; genes(g).exons{t}(1:end-1,2)+1, genes(g).exons{t}(2:end,1)-1];
  end
  introns = unique(introns, 'rows');
  genes(g).introns_covered = false(1,size(introns,1));
  genes(g).introns_conf = zeros(1,size(introns,1));
  try
    if isempty(intron_list) || isempty(introns)
      idx1 = []; idx2 = [];
    else
      [idx1 idx2] = ismember(introns,  intron_list(:,1:2), 'rows');
      idx1 = find(idx1);
      idx2 = idx2(idx1);
    end
  catch
    idx1 = []; idx2 = [];
  end
  genes(g).introns_covered(idx1) = true;
  genes(g).introns_conf(idx1) = intron_list(idx2,3);
end

[fd msg] = fopen(output_file, 'w+');
assert(fd~=-1);

%%%%% print result to stdout and output file %%%%%
used_contig_idx = unique([genes.chr_num]);
fprintf(fd, 'Created with ReadStats version %s\n', RQUANT_VERSION);
for f = [1 fd],
  % number of reads
  fprintf(f, '\n*** number of reads ***\n');
  for c = 1:length(contig_names),
    fprintf(f, 'contig %s: %i\n', contig_names{c}, mapped_reads(c));
  end
  fprintf(f, 'total: %i\n', sum(mapped_reads(used_contig_idx)));

  % median coverage
  mean_ec = [genes.mean_ec];
  fprintf(f, '\n*** Read coverage in %i genes ***\n', length(genes));
  fprintf(f, '25th percentile: %.2f \n', prctile(mean_ec, 25));
  fprintf(f, '50th percentile: %.2f \n', prctile(mean_ec, 50));
  fprintf(f, '75th percentile: %.2f \n', prctile(mean_ec, 75));

  % number of covered introns
  introns_covered = [genes.introns_covered];
  fprintf(f, '\n*** Number of covered introns ***\n');
  fprintf(f, '%i out of %i\n', sum(introns_covered), length(introns_covered));
  % median intron coverage
  introns_conf = [genes.introns_conf];
  fprintf(f, '\n*** Intron coverage of %i introns***\n', length(introns_covered)); 
  fprintf(f, '25th percentile: %.2f \n', prctile(introns_conf, 25));
  fprintf(f, '50th percentile: %.2f \n', prctile(introns_conf, 50));
  fprintf(f, '75th percentile: %.2f \n', prctile(introns_conf, 75));
end
  
fclose(fd);

[ret, timedate] = unix('date');
timedate(timedate==sprintf('\n')) = [];
fprintf(1, '\n-------------------------------------------------------------- \n');
fprintf(1, ' ReadStats version %s finished %s\n', RQUANT_VERSION, timedate) ; 
fprintf(1, '-------------------------------------------------------------- \n');
