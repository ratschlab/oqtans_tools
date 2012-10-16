function genes = build_splice_graph_caller(genes)
% BUILD_SPLICE_GRAPH_CALLER   Builds a splicegraph from for genes.
%
%   genes = build_splice_graph_caller(genes)
%
%   -- input --
%   genes: struct defining genes with start, stops, exons etc.
%
%   -- output --
%   genes: augmented by splicegraph
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2005-2009 Cheng Soon Ong, Gunnar Raetsch 
%   Copyright (C) 2005-2009 Max Planck Society
%

  
for gene_idx = 1:length(genes)
  % construct a new empty splice graph, based on
  % http://proline.bic.nus.edu.sg/dedb/
  vertices =  [];
  edges = [];
  if (mod(gene_idx,100)==0)
    fprintf(1,'.');
  end

  for transcript_idx = 1:length(genes(gene_idx).transcripts)
    exon_start_end = genes(gene_idx).exons{transcript_idx};
    if (size(exon_start_end,1) == 1)
      exon1_start = exon_start_end(1,1);
      exon1_end = exon_start_end(1,2);
      if isempty(vertices)
	vertices(1,1) = exon1_start;
	vertices(2,1) = exon1_end;
	edges = 0;
	num_exons = 1;
      else
	vertices(1,num_exons+1) = exon1_start;
	vertices(2,num_exons+1) = exon1_end;
	edges(1,num_exons+1) = 0;
	edges(num_exons+1,1) = 0;
	num_exons = num_exons + 1;
      end
    else
      for exon_idx = 1:(size(exon_start_end,1)-1)
	exon1_start = exon_start_end(exon_idx,1);
	exon1_end = exon_start_end(exon_idx,2);
	exon2_start = exon_start_end(exon_idx+1,1);
	exon2_end = exon_start_end(exon_idx+1,2);
	if isempty(vertices)
	  vertices(1,1) = exon1_start;
	  vertices(2,1) = exon1_end;
	  vertices(1,2) = exon2_start;
	  vertices(2,2) = exon2_end;
	  edges = zeros(2);
	  edges(1,2) = 1;
	  edges(2,1) = 1;
	  num_exons = 2;
	else
	  exon1_idx = 0;
	  exon2_idx = 0;
	  idx = 1;
	  while idx <= num_exons
	    if ((vertices(1,idx)==exon1_start) && (vertices(2,idx)==exon1_end))
	      exon1_idx = idx;
	    end
	    idx = idx+1;
	  end
	  idx = 1;
	  while idx <= num_exons
	    if ((vertices(1,idx)==exon2_start) && (vertices(2,idx)==exon2_end))
	      exon2_idx = idx;
	    end
	    idx = idx+1;
	  end
	  if (exon1_idx~=0) && (exon2_idx~=0)
	    edges(exon1_idx,exon2_idx) = 1;
	    edges(exon2_idx,exon1_idx) = 1;
	  else
	    if ((exon1_idx==0) && (exon2_idx~=0))
	      vertices(1,num_exons+1) = exon1_start;
	      vertices(2,num_exons+1) = exon1_end;
	      edges(exon2_idx,num_exons+1) = 1;
	      edges(num_exons+1,exon2_idx) = 1;
	      num_exons = num_exons + 1;
	    elseif ((exon2_idx==0) && (exon1_idx~=0))
	      vertices(1,num_exons+1) = exon2_start;
	      vertices(2,num_exons+1) = exon2_end;
	      edges(exon1_idx,num_exons+1) = 1;
	      edges(num_exons+1,exon1_idx) = 1;
	      num_exons = num_exons + 1;
	    else
	      assert((exon1_idx==0)&&(exon2_idx==0));
	      vertices(1,num_exons+1) = exon1_start;
	      vertices(2,num_exons+1) = exon1_end;
	      num_exons = num_exons + 1;
	      vertices(1,num_exons+1) = exon2_start;
	      vertices(2,num_exons+1) = exon2_end;
	      num_exons = num_exons + 1;	  
	      
	      edges(num_exons-1,num_exons) = 1;
	      edges(num_exons,num_exons-1) = 1;
	    end
	  end
	end
      end
    end
  end
  genes(gene_idx).splicegraph = {vertices,edges};
end
fprintf(1,'\n');