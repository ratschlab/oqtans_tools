function mTIM_galaxy_train(...
    genome_dir, acc_splice_dir, don_splice_dir, read_map_file, anno_dir, out_dir)
% 
% Galaxy configuration script.
%
% Inputs:
%   genome_dir        : path and filename of .gio file 
%   acc_splice_dir    : path to the acceptor splice site predictions  
%   don_splice_dir    : path to donor splice site predictions
%   read_map_file     : path and filename of .bam read-maps
%   anno_dir          : path to the annotation (a genes.mat file
%                       is assumed in this directory) 
%   out_dir           : the directory where all the output is written to
%                       (data preparation files, training files, predictions, ...)
%
%
% written by Nico Goernitz, Berlin Institute of Technology, 2012

% setup the cofiguration
CFG = general_settings(out_dir,'training',1,1);

% complete the paths
% genome details
CFG.genome_dir = genome_dir;
CFG.genome_info = sprintf('%s/genome.config', genome_dir);
info = init_genome(CFG.genome_info);
CFG.chr_names = info.contig_names;
CFG.num_chr = length(CFG.chr_names);
for c=1:CFG.num_chr,
  info.flat_fnames{c} = [info.basedir '/genome/' info.contig_names{c} '.flat'];
  d = dir(info.flat_fnames{c});
  CFG.chr_lens(c) = d.bytes;
  assert(CFG.chr_lens(c)>0);
end

% annotation details
CFG.annotation_dir = sprintf('%s', anno_dir);
CFG.gene_fn = sprintf('%s/genes.mat', anno_dir);
% splice site prediction details
CFG.splice_site_dir.acc = acc_splice_dir;
CFG.splice_site_dir.don = don_splice_dir;
% get the read map directory
idx = find(read_map_file=='/', 1, 'last');
CFG.read_map_dir = read_map_file(1:idx);
CFG.read_map_file = read_map_file;
fprintf('dir:%s\nfile:%s\n',CFG.read_map_dir,CFG.read_map_file);

% 1. step 
% generate data
% - divide into train/vald/test data chunks
% - fill chunks with features
main_generate_data(CFG);
fprintf('Data Generated.\n');

% 2. step
% training
% - train a classifier using the prepared training
%   and validation chunks
main_run_training(CFG);
fprintf('Training finished.\n');
