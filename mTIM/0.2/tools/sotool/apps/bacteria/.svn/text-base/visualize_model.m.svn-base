function visualize_model(modelPath)
% 
%
%
clc

path = 'model';
addpath(path);

% get function names
config = model_config();

% get the state model
PAR = [];
PAR.num_features = 1;

stateModel = eval([config.func_make_model '(PAR)']);
% number of states
states = length(stateModel);
rmpath(path);

text = sprintf('digraph G {\n');

% principles
for i=1:length(stateModel),
  name = stateModel(i).name;

  id = stateModel(i).id;
  name = stateModel(i).name;
  is_start = stateModel(i).is_start;
  is_stop = stateModel(i).is_stop;
  label = stateModel(i).label;

  successors = stateModel(i).successors;          % e.g. [3 4]
  trans_scores = stateModel(i).trans_scores;      % e.g. [0 0]
  learn_scores = stateModel(i).learn_scores;      % [4x1 double]
  feature_scores = stateModel(i).feature_scores;  % [0x2 double]
  monot_scores = stateModel(i).monot_scores;      % [0x1 double]
  score_coupling = stateModel(i).score_coupling;  % [0x2 double]  
  
  if (stateModel(i).is_start),
    text = [text sprintf('%s %s;\n',name,'[shape=box,color=green]')];
  elseif (stateModel(i).is_stop),
    text = [text sprintf('%s %s;\n',name,'[shape=box,color=red]')];
  else
    text = [text sprintf('%s %s;\n',name,'[shape=ellipse]')];
  end
end  
  
for i=1:length(stateModel),
  node = stateModel(i);
  name1 = node.name;
  chid = node.successors(find(node.successors~=0));
  
  successors = stateModel(i).successors;          % e.g. [3 4]
  trans_scores = stateModel(i).trans_scores;      % e.g. [0 0]

  learn_scores = stateModel(i).learn_scores;      % [4x1 double]
  feature_scores = stateModel(i).feature_scores;  % [0x2 double]
  monot_scores = stateModel(i).monot_scores;      % [0x1 double]
  score_coupling = stateModel(i).score_coupling;  % [0x2 double]  

  for k=1:length(chid),
    j = chid(k);
    name2 = stateModel(j).name;
    
    add = '[';
    if (trans_scores(k)>0), add = [add 'style=bold'];
    else add = [add 'color=lightgrey']; end
    add = [add ']'];
    text = [text sprintf('%s->%s %s;\n',name1,name2,add)];
  end
end
text = [text '}'];
fid = fopen('graph.dot','w');
fwrite(fid,text);
fclose(fid);

text
system('dot -Tpdf graph.dot -o test.pdf');
system('evince test.pdf');


