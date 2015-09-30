function [models,param_names] = generate_all_models(Model_param)
% models = generate_all_models(Model_param) 
  
param_names = fieldnames(Model_param); 
if isempty(param_names)
  models =[];param_names=[];
  return
end

for j = 1:length(param_names)
  param_nums(j) = length(getfield(Model_param,param_names{j}));
end
models = zeros(length(param_names),prod(param_nums));
for j = 1:length(param_names)
  par_tmp  = getfield(Model_param,param_names{j});
  par1 = [];
  for k = 1:length(par_tmp)
    par1 = [par1 repmat(k,1,prod(param_nums(j+1:end)))];
  end
  models(j,:) = par_tmp(repmat(par1,1,prod(param_nums(1:j-1))));
end

if prod(param_nums)<1
  assert(length(unique(models','rows'))==prod(param_nums))
else
  assert(~isempty(models))
end
