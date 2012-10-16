function plot_evaluation()

% TODO hacky
addpath ../utils

methods = {'mtim', 'cufflinks'};%, 'scripture'};

%data_dir = '/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/mTIM/';
data_dir = '/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana/';

%exp = 'elegans.6';
%exp = 'elegans.6.bestpaired.filtered';
%exp = 'elegans.tophat.PE';
%exp = 'elegans.tophat.PE.bestpaired.filtered';
exp = 'SRX_all.new'

for m=1:length(methods),
  eval_dir = sprintf('%smtim_eval/%s/%s/', data_dir, exp, methods{m});
  try
    fn = strtrim(ls(eval_dir));
  catch
    warning('unable to load evaluation information for %s from %s', ...
            methods{m}, eval_dir);
    continue
  end
  % TODO assert that this is only a single file!
  fns = splitstr(fn);
  assert(length(fns) == 1);
  load([eval_dir fns{1}], 'eval');
  eval.method = methods{m};
  methods_eval(m) = eval;
end

bh = zeros(3, length(methods_eval));
for m=1:length(methods_eval),
  lg{m} = methods_eval(m).method;
  bh(1,m) = methods_eval(m).intron.f1;
  bh(2,m) = methods_eval(m).transcr.f1;
  bh(3,m) = methods_eval(m).gene.f1;
end

figure();
bar(bh, 'grouped');
title(sprintf('Accuracy on %s', exp));
ylim([0, 1])
ylabel('F1-Score');
set(gca, 'XTick', []);
th = text(1, 0.9, 'intron');
set(th, 'HorizontalAlignment', 'center');
th = text(2, 0.9, 'transcript');
set(th, 'HorizontalAlignment', 'center');
th = text(3, 0.9, 'gene');
set(th, 'HorizontalAlignment', 'center');
legend(lg, 'Location', 'East');

fig_dir =  sprintf('%smtim_eval/%s/', data_dir, exp);
dt = datestr(now, 'yyyy-mm-dd_HHhMM');
print('-depsc', [fig_dir 'comparative_eval_' dt '.eps'])

keyboard

