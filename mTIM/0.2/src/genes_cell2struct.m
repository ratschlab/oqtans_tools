function genes_cell2struct(out_dir)
% GENES_CELL2STRUCT   Converts genes stored as a cell to struct.
%
%   genes_cell2struct(anno_fname)
%
%   -- input --
%   anno_fname:   name of file where genes as cell are stored
%
%   -- output --
%   genes as a struct

anno_fname = strcat(out_dir, '/genes.mat');
genome_info_file = strcat(out_dir, '/genome.config');
genome_info = init_genome(genome_info_file);
load(anno_fname, 'genes');
if iscell(genes)
  genes_cell = genes;
  clear genes;
  for g = 1:length(genes_cell), 
    gene = genes_cell{g};
    for e = 1:length(gene.exons)
      gene.exons{e} = double(gene.exons{e});
    end    
    gene.exons = reshape(gene.exons, 1, length(gene.exons));
    gene.id = double(gene.id);
    gene.start = double(gene.start);
    gene.stop = double(gene.stop);
    gene.chr_num = strmatch(upper(gene.chr), upper(genome_info.contig_names), 'exact');
    genes(g) = gene;
  end
  save(anno_fname, 'genes');
end
