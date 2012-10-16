function write_rquant_gff3(CFG, genes, source, gff_fname, mapped_reads) 
% WRITE_RQUANT_GFF3   Writes rQuant results in gff3 format.
%
% write_rquant_gff3(CFG, genes, source, gff_fname, mapped_reads)
%
%   -- input --
%   CFG:          configuration struct
%   genes:        struct defining genes with start, stops, exons etc.
%   source:       tool with which the data was generated
%   gff_fname:    name of gff file
%   mapped_reads: number of mapped reads
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


used_contig_idx = unique([genes.chr_num]);
mapped_reads = sum(mapped_reads(used_contig_idx));
reads_per_1coverage = 1000/CFG.read_len;
rpkm_factor_all = reads_per_1coverage./(mapped_reads/10^6);
rpkm_factor_weighting = mapped_reads.*CFG.read_len;
rpkm_factor_weighting = rpkm_factor_weighting/sum(rpkm_factor_weighting);
rpkm_factor = sum(rpkm_factor_all.*rpkm_factor_weighting);

if CFG.VERBOSE>1
  if exist(gff_fname, 'file')
    fprintf(1, 'replacing file %s...\n', gff_fname);
  else
    fprintf(1, 'creating file %s...\n', gff_fname);
  end
end
[fd msg] = fopen(gff_fname, 'w+');
if CFG.VERBOSE>1, disp(msg); end
assert(fd~=-1);

fprintf(fd, '##gff-version 3\n');
for g = 1:length(genes),
  if CFG.VERBOSE>1, fprintf(1, 'writing gene %i...\r', g); end
  gene = genes(g);
  type  = 'gene';
  score = '.';
  phase = '.';
  start = gene.start;
  stop = gene.stop;
  attr_str = sprintf('ID=%s;Name=%s', gene.name, gene.name);
  fprintf(fd, '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n', gene.chr, source, type, start, stop, score, gene.strand, phase, attr_str);
  for t = 1:length(gene.transcripts),
    type = 'transcript';
    score = '.';
    phase = '.';
    start = min(gene.exons{t}(:,1));
    stop = max(gene.exons{t}(:,2));
    % transcript
    if isfield(gene, 'transcript_weights'), % write RPKM to attribute field
      attr_str = sprintf('ID=%s;Parent=%s;ARC=%4.4f;RPKM=%4.4f', gene.transcripts{t}, gene.name, gene.transcript_weights(t), gene.transcript_weights(t)*rpkm_factor);
    else
      attr_str = sprintf('ID=%s;Parent=%s', gene.transcripts{t}, gene.name);
    end
    fprintf(fd, '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n', gene.chr, source, type, start, stop, score, gene.strand, phase, attr_str);
    % exons
    attr_str = sprintf('Parent=%s', gene.transcripts{t});
    exon_type = {'exons'};
    gff_types  = {'exon'};
    for tt = 1:length(exon_type),
      exons = gene.(exon_type{tt}){t};
      for e = 1:size(exons,1),
        type = gff_types{tt};
        score = '.';
        phase='.';
        start = exons(e,1);
        stop = exons(e,2);
        fprintf(fd, '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n', gene.chr, source, type, start, stop, score, gene.strand, phase, attr_str);
      end
    end
  end
end
fclose(fd);
