function feats = append_repeat_feats(chunks, exp, feats, CFG)

% feats = append_repeat_feats(chunks, exp, feats, CFG)
%
% 
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009


% TODO this function is probably outdated/nonfunctional!
warning('Calling outdated/nonfucntional append_repeat_feats.m function!')


assert(issorted(chunks, 'rows'));

strands = '+-';
% retrieve position and repeat flags by interval_query 
% (see interval_query.c in ~/svn/projects/genefinding/utils)
fields = {'repeats'};

flag_map = [32 96 224; 1 2 3];

L = sum(chunks(:,3)-chunks(:,2)+1);
repeat_feat = zeros(1, L);

data_dir = sprintf('%s%s_repeat_maps/%s_%s_repeat', CFG.repeat_dir, exp);
if ~exist(data_dir, 'dir'),
  fprintf('  no repeat information found for %s, skipping ...\n', exp);
  feats = [feats; repeat_feat];
  return
end
fprintf('  collecting repeat data for %s ...\n', exp);

offset = 0;
%tic
for c=1:size(chunks,1),
  chr = chunks(c,1);
  % collect repeat data
  fn = sprintf('%s%s_%s_repeat', data_dir, exp, CFG.chr_names{chr});
  if exist([fn '.pos'], 'file') ~= 0,
    [pos, flag] = interval_query(fn, fields, [chunks(c,2); chunks(c,3)]);
    assert(all(chunks(c,2)<=pos & pos<=chunks(c,3)));
    pos = pos' - chunks(c,2) + 1 + offset;
    for f=1:size(flag_map,2),
      idx = flag==flag_map(1,f);
      repeat_feat(pos(idx)) = flag_map(2,f);
    end
  else
%    warning('Could not find %s! (ignoring)', [fn '.pos'])
  end  
  offset = offset + chunks(c,3)-chunks(c,2)+1;
%  if mod(c,100)==0,
%    fprintf('    processed %i chunks (%2.1f%%, %.1f sec)\r', ...
%            c, 100*c/size(chunks,1), toc);
%  end
end
%fprintf('    processed %i chunks (%2.1f%%, %.1f sec)\n', ...
%        c, 100*c/size(chunks,1), toc);

feats = [feats; repeat_feat];
assert(size(feats,2) == L);

% eof