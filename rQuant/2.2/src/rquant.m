function rquant(anno_dir, track, output_file, output_dir, load_profiles, profiles_fn, learn_profiles, profiles_fn_out, CFG)
% RQUANT   Determines the abundance of multiple transcripts per gene locus from RNA-Seq measurements.
%
%   rquant(anno_dir, track, output_file, output_dir, load_profiles, profiles_fn, learn_profiles, profiles_fn_out, CFG)
%
%   -- input --
%   anno_dir:        directory of genes
%   track:           name of BAM file
%   output_file:     result gff3 file 
%   output_dir:      output directory
%   load_profiles:   flag that enables loading of profiles
%   profiles_fn:     name of input file that stores profiles
%   learn_profiles:  flag that enables learning of profiles
%   profiles_fn_out: name of output file that stores profiles 
%   CFG:             configuration struct (optional, allows to pass parameters settings)
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


if exist('load_profiles', 'var')
  if ~isnumeric(load_profiles)
    load_profiles = str2num(load_profiles);
  end
end
if exist('learn_profiles', 'var')
  if ~isnumeric(learn_profiles)
    learn_profiles = str2num(learn_profiles);
  end
end


%%%%% options Galaxy/local usage %%%%%
% more output to stdout
CFG.VERBOSE = 1; % 0: no output, 1: more output, 2: debug output
% parallelisation with rproc (1: cluster submission or 0: locally)
if ~isfield(CFG, 'use_rproc'), CFG.use_rproc = 0; end


%%%%% paths %%%%%
% rQuant paths
global RQUANT_PATH RQUANT_SRC_PATH
% rQuant version
global RQUANT_VERSION
% SAMTools path
global SAMTOOLS_DIR
CFG.paths = sprintf('%s/mex:%s/tools:%s', RQUANT_PATH, RQUANT_PATH, RQUANT_SRC_PATH);
if CFG.use_rproc
  CFG.paths = sprintf('%s:~/svn/tools/rproc', CFG.paths);
end
addpath(CFG.paths);


%%%%% output directory %%%%%
CFG.out_dir = sprintf('%s/', output_dir);
if ~exist(CFG.out_dir ,'dir'),
  [s m mid] = mkdir(CFG.out_dir);
  assert(s);
end
CFG.output_file = strrep(output_file, '.gff3', '.mat');


%%%%% directories from which to load read data and genes %%%%%
CFG.samtools_dir = sprintf('%s/', SAMTOOLS_DIR);
CFG.gene_dir = anno_dir; 
CFG.gene_fn = sprintf('%s/genes.mat', CFG.gene_dir);
if ~isfield(CFG, 'repeats_fn'), CFG.repeats_fn = ''; end


%%%%% genome information %%%%%
try
  if iscell(track)
    [s bam_header] = unix(sprintf('%s./samtools view %s -H', CFG.samtools_dir, track{1}));
  else
    [s bam_header] = unix(sprintf('%s./samtools view %s -H', CFG.samtools_dir, track));
  end
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
if ~isfield(CFG, 'read_len')
  assert(length(track) & ~iscell(track));
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
else
  for c = 1:length(contig_names),
    if iscell(track)
      for f = 1:length(track),
        CFG.tracks_fn{c}{f} = track{f};
      end
    else
      CFG.tracks_fn{c} = {track};
    end
  end
  mapped_reads = zeros(1, length(CFG.tracks_fn));
end


%%%%% rquant parameters %%%%%

% enables taking data for both strands together
CFG.both_strands = 1;

%%%%% transcript weight optimisation
% method to determine transcript weights 
CFG.method = 'pos'; % 'pos' or 'seg'
if ~isfield(CFG, 'correct_intervals'), CFG.correct_intervals = 0; end % correction to closed interval
if ~isfield(CFG, 'C_I'), CFG.C_I = 100; end  
if ~isfield(CFG, 'paired'), CFG.paired = 0; end % usage of paired-end data
if CFG.paired && ~isfield(CFG, 'C_PE'), CFG.C_PE = 100; end
if CFG.paired, CFG.ins_sizes = []; end

%%%%% profile learning
% enables loading of profiles from CFG.profiles_fn
CFG.load_profiles = load_profiles;
CFG.profiles_fn = profiles_fn;
CFG.learn_profiles = learn_profiles; % 0: no learning, 1: empirically estimated, 2: optimised
% number of iterations
CFG.max_iter = 100;
% number of supporting points for profile functions
CFG.num_plifs = 50;
% maximal number of positions to be considered at both transcript ends
CFG.max_side_len = 500;
% bins for different transcript lengths
CFG.transcript_len_ranges = [];
% enables subsampling of data for learning profiles
CFG.subsample = 0;
% maximal number of examples for learning profiles
CFG.max_num_train_exm = 1e6;
% fraction of positions to be subsampled for learning profiles
CFG.subsample_frac = 0.2;
% regularisation strengths
if ~isfield(CFG, 'C_F'), CFG.C_F = 10^4; end
if ~isfield(CFG, 'C_N'), CFG.C_N = 10; end
% checks of variable domains
if CFG.load_profiles && CFG.learn_profiles==1,
  error('Pre-learned profiles cannot be used for empirical profile estimation.');
end
if ~exist(profiles_fn, 'file') && CFG.load_profiles
  error('File with pre-learned profiles does not exist.');
end


%%%%% rquant %%%%%
save_fname = rquant_core(CFG);


%%%%% write to gff file %%%%%
if ~isfield(CFG, 'write_gff'), CFG.write_gff = 1; end
if CFG.write_gff && ~isempty(output_file)
  load(save_fname, 'genes');
  write_rquant_gff3(CFG, genes, sprintf('rQuant v%s', RQUANT_VERSION), output_file, mapped_reads);
end

%%%%% write learned read density model %%%%%
if ~isfield(CFG, 'write_density_model'), CFG.write_density_model = 1; end
if CFG.write_density_model && CFG.learn_profiles>0 && ~isempty(profiles_fn_out)
  load(save_fname, 'profile_weights');
  write_density_model(CFG, profile_weights, profiles_fn_out);
end
