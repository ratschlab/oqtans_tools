function models= append_new_models(new_models, old_models)

% models= append_new_models(new_models, old_models)
% 
models = old_models;
 
for j=1:size(new_models,2)
  model = new_models(:,j);
  idx=1:size(old_models,2);
  for k=1:length(model)
    idx = intersect(idx, find(old_models(k,:)==model(k)));
  end
  if isempty(idx)
    models = [models model];
  else
    assert(length(idx)==1)
  end
end
