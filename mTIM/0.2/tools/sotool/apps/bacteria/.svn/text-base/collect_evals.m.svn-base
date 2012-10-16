function collect_evals(ev_files)
%
%
% adapted by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011

fprintf('\nEvaluate %i files.\n',length(ev_files));
fprintf('Only print results of files that are either leafs or root nodes.\n');

for i=1:length(ev_files),
  % partition data for cross-validation
  data_file = [ev_files{i}];
  load(data_file);

  % flag for root nodes
  isRoot = ~PAR.mtl.mtl_enable && length(PAR.organism)>1;
  % flag for leafs
  isLeaf = PAR.mtl.mtl_enable && length(PAR.organism)==1;
  
  if (isLeaf),
    mtlfile = PAR.mtl.filename;
    fprintf('Leaf %s(#train=%i): MTL enabled with file (%s) b,C=%1.2f,%3.2f.\n',data_file,length(PAR.train_exms),mtlfile,PAR.mtl.mtl_b,PAR.C_small);
    % only one organism
    fprintf('Organism(#train): %s(%i) with f_score %1.4f.\n',res{1}.name,length(PAR.organism{1}.train_exms),res{1}.result.f_score);
    fprintf('SINGLE Organism(#train): %s(%i) with f_score %1.4f.\n',res_single.name,length(PAR.organism{1}.train_exms),res_single.result.f_score);
  end
  
  if (isRoot),
    fprintf('Root-node %s(#train=%i): MTL disabled C=%3.2f.\n',data_file,length(PAR.train_exms),PAR.C_small);
    % print results of all organisms
    for j=1:length(res),
      fprintf('Organism(#train): %s(%i) with f_score %1.4f.\n',res{j}.name,length(PAR.organism{j}.train_exms),res{j}.result.f_score);
    end
  end
    
  fprintf('\n');
end
% eof

