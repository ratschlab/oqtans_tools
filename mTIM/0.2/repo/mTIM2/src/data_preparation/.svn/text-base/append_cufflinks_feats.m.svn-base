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
    warning('Cufflinks feature will be empty (noo file set)!');
    feats = [feats; zeros(1,L)];
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
cf = zeros(1,L);


for c=1:size(chunks,1),
  region.start   = chunks(c,2);
  region.stop    = chunks(c,3);
  region.chr_num = chunks(c,1);
  region.chr = CFG.chr_names{chunks(c,1)};
  region.len = region.stop-region.start+1;
  
  % find all genes within the current region
  inds = find(strcmpi({cuffl.genes.chr},region.chr));
  % build exon list
  exon_list = [];
  for i=1:length(inds),
      for j=1:length(cuffl.genes(inds(i)).exons),
          exons = cuffl.genes(inds(i)).exons{j};          
          for k=1:size(exons,1),
            % TODO: CHECK!
            start = exons(k,1) - region.start;
            stop  = exons(k,2) - region.start;
            
            % only handle cases where start>0 and stop<len
            if (start>0 && stop<region.len),
                idxs = offset+start+1 : offset+stop+1;
                cf(1,idxs) = cf(1,idxs) + 1;
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
