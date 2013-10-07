function feats = append_cufflinks_feats(chunks, feats, CFG)
% feats = append_cufflinks_feats(chunks, feats, CFG)
%
% Most basic version:
% - cufflinks coverage per positions (integers)
% - one feature for both strands
%
% written by Nico Görnitz, 04-2012, Berlin Institute of Technology

assert(issorted(chunks, 'rows'));
tic

L = sum(chunks(:,3)-chunks(:,2)+1);
% features correspond to a horizontal concatenation (of length L) of feature
% blocks of all individual chunks

strands = '+-';
offset = 0;

% check if cufflinks prediction file is set
if (isfield(CFG,'cufflinks_pred_file') && exist(CFG.cufflinks_pred_file,'file')),
    fprintf('  using cufflinks prediction file as feature.');
else
    warning('Cufflinks features will be empty (noo file set)!');
    feats = [feats; zeros(2,L)];
    return;
end

% load cufflinks predicted genes structure
cuffl = [];
fprintf('Loading cufflinks genes from "%s"...\n', CFG.cufflinks_pred_file);
cuffl = load(CFG.cufflinks_pred_file, 'genes');
fprintf('  loaded %i predicted genes.\n', length(cuffl.genes));
fprintf('  converting genes: %i\n',CFG.cufflinks_convert);
if (CFG.cufflinks_convert),
    cuffl.genes = closed_to_half_open(cuffl.genes);
end



% feature sequence
% 0 : IGE 
% 1 : Intron 
% 2 : Exon  
% strand W=1 C=2
cf = zeros(2,L);


for c=1:size(chunks,1),
  region.start   = chunks(c,2);
  region.stop    = chunks(c,3);
  region.chr_num = chunks(c,1);
  region.chr = CFG.chr_names{chunks(c,1)};
  region.len = region.stop-region.start+1;
  
  % find all genes within the current region
  inds = find(strcmpi({cuffl.genes.chr},region.chr));
  % genes that lie inside the area
  n1inds = find([cuffl.genes(inds).start]>=region.start & [cuffl.genes(inds).stop]<region.stop);
  % genes that end in the area
  n2inds = find([cuffl.genes(inds).start]<region.start & [cuffl.genes(inds).stop]<region.stop);
  % genes that start in the area
  n3inds = find([cuffl.genes(inds).start]>=region.start & [cuffl.genes(inds).stop]>region.stop);
  % and genes that completely cover the area
  n4inds = find([cuffl.genes(inds).start]<region.start & [cuffl.genes(inds).stop]>region.stop);

  inds = inds([n1inds, n2inds, n3inds, n4inds]);

  % build exon list
  exon_list = [];
  for i=1:length(inds),

      strand = cuffl.genes(inds(i)).strand;
      strand_ind = 1;
      if strcmpi(strand,'-'), strand_ind=2; end;

      % mark the whole gene as either intron W or intron C
      start = cuffl.genes(inds(i)).start - region.start;
      if (start<0), start=0; end;
      stop = cuffl.genes(inds(i)).stop - region.start;
      if (stop>=region.len), stop=region.len-1; end;
      idxs = offset+start+1 : offset+stop+1;
      cf(1,idxs) = 1;

      % overwrite exonic regions with exon W or exon C
      for j=1:length(cuffl.genes(inds(i)).exons),
          exons = cuffl.genes(inds(i)).exons{j};          
          for k=1:size(exons,1),
            % TODO: CHECK!
            start = exons(k,1) - region.start;

            if (start<0), start=0; end;
            if (start>=region.len), continue; end;

            stop  = exons(k,2) - region.start;
            if (stop<0), continue; end;
            if (stop>=region.len), stop=region.len-1; end;
            
            % only handle cases where start>0 and stop<len
            if (start>0 && stop<region.len),
                idxs = offset+start+1 : offset+stop+1;
                cf(strand_ind,idxs) = 2;
            end
          end
      end
  end

  % advance offset to point to the start of the next feature chunk
  offset = offset + region.stop - region.start + 1;
end


feats = [feats; cf];
assert(size(feats,2) == L);

fprintf('  appended cufflinks features to %i chunks in %.1f sec\n', size(chunks,1), toc);




% eof
