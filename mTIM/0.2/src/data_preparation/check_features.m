function check_features(signal, label, sample_size)

% check_features(signal, label, sample_size)
%
% Checks the predictivity of features using ROC and PRC analysis.
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009-2011

LABELS = get_label_set_mTIM();

fprintf(['Composition of labelings in %12i genomic positions included in ' ...
         'training set \n'], length(label));
fn = fieldnames(LABELS);
for f=1:length(fn),
  cnt = sum(label == getfield(LABELS, fn{f}));
  if cnt > 0,
    fprintf('Labeled %12i training positions as %10s (%2.1f%%)\n', ...
            cnt, fn{f}, 100*cnt/length(label));
  end
end

% sample positions
s = 1:length(label);
s(label == LABELS.ambiguous) = [];
s = s(randperm(length(s)));
s = s(1:sample_size);

gidx = label==LABELS.exon_W | label==LABELS.exon_C ...
       | label==LABELS.intron_W | label==LABELS.intron_C;

%keyboard

% exon coverage feature predictivity
F = 1;
l = -ones(size(s));
l(label(s)==LABELS.exon_W | label(s)==LABELS.exon_C) = +1;
x = symlog2(signal(F,s));
[TP FP TN FN] = eval_separation(x, l, 50);
OPT.line_style = {'.'};
h = figure
plot_PRC(TP, FP, TN, FN, OPT);
title('exon coverage feature predictivity, PRC');
grid on;
saveas(h,'exon_coverage.eps','eps');
% exon difference feature predictivity
F = 2;
l = -ones(size(s));
l(label(s)==LABELS.exon_W) = +1;
idx = find(label(s)==LABELS.exon_W | label(s)==LABELS.exon_C);
l = l(idx);
x = symlog2(signal(F,s(idx)));
[TP FP TN FN] = eval_separation(x, l, 50);
OPT.line_style = {'.'};
h = figure
plot_PRC(TP, FP, TN, FN, OPT);
title('exon difference feature predictivity, PRC');
grid on;
saveas(h,'exon_difference.eps','eps');

% exon coverage gradient feature predictivity
F = 3;
ebw = find_blocks(label == LABELS.exon_W);
ebc = find_blocks(label == LABELS.exon_C);
l = -ones(size(label));
l(ebw(1,:)) = +1;
l(ebc(1,:)) = +1;
l(ebw(2,:)) = +1;
l(ebc(2,:)) = +1;
lambda = 2;
for i=1:lambda,
  l(ebw(1,:)-i) = +1;
  l(ebw(1,:)+i) = +1;
  l(ebw(2,:)-i) = +1;
  l(ebw(2,:)+i) = +1;
  l(ebc(1,:)-i) = +1;
  l(ebc(1,:)+i) = +1;
  l(ebc(2,:)-i) = +1;
  l(ebc(2,:)+i) = +1;
end
x = signal(F,:);
x = abs(symlog2(x));
fprintf('eval_separation..');
[TP FP TN FN] = eval_separation(x, l, 50);
fprintf('finished.\n');
OPT.line_style = {'.'};
h = figure
plot_PRC(TP, FP, TN, FN, OPT);
title('exon coverage gradient feature predictivity wrt. exon boundaries, PRC');
grid on;
saveas(h,'exon_coverage_gradient.eps','eps');


% intron span feature predictivity
F = 4;
l = -ones(size(s));
l(label(s)==LABELS.intron_W | label(s)==LABELS.intron_C) = +1;
x = symlog2(signal(F,s));
[TP FP TN FN] = eval_separation(x, l, 50);
OPT.line_style = {'.'};
h = figure
plot_PRC(TP, FP, TN, FN, OPT);
title('intron span feature predictivity, PRC');
grid on;
saveas(h,'intron_span.eps','eps');

% intron difference feature predictivity
F = 5;
l = -ones(size(s));
l(label(s)==LABELS.intron_W) = +1;
idx = find(label(s)==LABELS.intron_W | label(s)==LABELS.intron_C);
l = l(idx);
x = symlog2(signal(F,s(idx)));
[TP FP TN FN] = eval_separation(x, l, 50);
OPT.line_style = {'.'};
h=figure
plot_PRC(TP, FP, TN, FN, OPT);
title('intron difference feature predictivity, PRC');
grid on;
saveas(h,'intron_difference.eps','eps');

ibw = find_blocks(label == LABELS.intron_W);
ibc = find_blocks(label == LABELS.intron_C);

% read intron start feature W strand
F = 6;
l = -ones(size(label));
l(ibw(1,:)) = +1;
x = signal(F,:);
x(isnan(x)) = 0;
[TP FP TN FN] = eval_separation(x, l);
OPT.line_style = {'.'};
h=figure
plot_PRC(TP, FP, TN, FN, OPT);
title('intron start feature predictivity (+ strand), PRC');
grid on;
saveas(h,'read_intron_start.eps','eps');

% read intron end feature W strand
F = 7;
l = -ones(size(label));
l(ibw(2,:)) = +1;
x = signal(F,:);
[TP FP TN FN] = eval_separation(x, l);
OPT.line_style = {'.'};
h=figure
plot_PRC(TP, FP, TN, FN, OPT);
title('intron end feature predictivity (+ strand), PRC');
grid on;
saveas(h,'read_intron_end.eps','eps');

% read intron start feature C strand
F = 8;
l = -ones(size(label));
l(ibc(1,:)) = +1;
x = signal(F,:);
[TP FP TN FN] = eval_separation(x, l);
OPT.line_style = {'.'};
h=figure
plot_PRC(TP, FP, TN, FN, OPT);
title('intron start feature predictivity (- strand), PRC');
grid on;
saveas(h,'read_intron_start_C.eps','eps');

% read intron end feature C strand
F = 9;
l = -ones(size(label));
l(ibc(2,:)) = +1;
x = signal(F,:);
[TP FP TN FN] = eval_separation(x, l);
OPT.line_style = {'.'};
h=figure
plot_PRC(TP, FP, TN, FN, OPT);
title('intron end feature predictivity (- strand), PRC');
grid on;
saveas(h,'read_intron_end_C.eps','eps');

% donor prediction feature W strand
F = 10;
l = -ones(size(label));
l(ibw(1,:)) = +1;
l = l(s);
x = signal(10,s);
x(x < 0) = 0;
[TP FP TN FN] = eval_separation(x, l, 250);
OPT.line_style = {'.'};
h=figure
plot_PRC(TP, FP, TN, FN, OPT);
title('DON feature predictivity (+ strand), PRC');
grid on;
saveas(h,'donor_pred_W.eps','eps');


% acceptor prediction feature W strand
F = 11;
l = -ones(size(label));
l(ibw(2,:)) = +1;
l = l(s);
x = signal(F,s);
x(x < 0) = 0;
[TP FP TN FN] = eval_separation(x, l, 250);
OPT.line_style = {'.'};
h=figure
plot_PRC(TP, FP, TN, FN, OPT);
title('ACC feature predictivity (+ strand), PRC');
grid on;
saveas(h,'acc_pred_W.eps','eps');

% donor prediction feature C strand
F = 12;
l = -ones(size(label));
l(ibc(2,:)) = +1;
l = l(s);
x = signal(F,s);
x(x < 0) = 0;
[TP FP TN FN] = eval_separation(x, l, 250);
OPT.line_style = {'.'};
h=figure
plot_PRC(TP, FP, TN, FN, OPT);
title('DON feature predictivity (- strand), PRC');
grid on;
saveas(h,'don_pred_C.eps','eps');

% acceptor prediction feature C strand
F = 13;
l = -ones(size(label));
l(ibc(1,:)) = +1;
l = l(s);
x = signal(F,s);
x(x < 0) = 0;
[TP FP TN FN] = eval_separation(x, l, 250);
OPT.line_style = {'.'};
h=figure
plot_PRC(TP, FP, TN, FN, OPT);
title('ACC feature predictivity (- strand), PRC');
grid on;
saveas(h,'acc_pred_C.eps','eps');


if size(signal,1) > F,
  % mate pair feature predictivity
  F = 14;
  l = ones(size(s));
  l(label(s)==LABELS.intergenic) = -1;
  x = symlog2(signal(F,s));
  [TP FP TN FN] = eval_separation(x, l, 50);
  OPT.line_style = {'.'};
 h= figure
  plot_ROC(TP, FP, TN, FN, OPT);
  title('mate pair feature predictivity wrt. exons & introns, ROC');
  grid on;
saveas(h,'mate_pair.eps','eps');
  h=figure
  plot_PRC(TP, FP, TN, FN, OPT);
  title('mate pair feature predictivity wrt. exons & introns, PRC');
  grid on;
  
saveas(h,'mate_pair_prc.eps','eps');
  % total coverage feature predictivity
  F = 15;
  l = ones(size(s));
  l(label(s)==LABELS.intergenic) = -1;
  x = symlog2(signal(F,s));
  [TP FP TN FN] = eval_separation(x, l, 50);
  OPT.line_style = {'.'};
  h=figure
  plot_ROC(TP, FP, TN, FN, OPT);
  title('total coverage feature predictivity wrt. exons & introns, ROC');
  grid on;
saveas(h,'total_cov.eps','eps');
  h=figure
  plot_PRC(TP, FP, TN, FN, OPT);
  title('total coverage feature predictivity wrt. exons & introns, PRC');
  grid on;
saveas(h,'total_cov_prc.eps','eps');

  % low-coverage block length feature predictivity
  F = 16;
  l = -ones(size(s));
  l(label(s)==LABELS.intergenic) = +1;
  x = signal(F,s);
  [TP FP TN FN] = eval_separation(x, l, 50);
  OPT.line_style = {'.'};
  h=figure
  plot_ROC(TP, FP, TN, FN, OPT);
  title('zero-cover block feature predictivity wrt. intergenic regions, ROC');
  grid on;
saveas(h,'zero_cov.eps','eps');
  h=figure
  plot_PRC(TP, FP, TN, FN, OPT);
  title('zero-cover block feature predictivity wrt. intergenic regions, PRC');
  grid on;
saveas(h,'zero_cov_prc.eps','eps');
end

if size(signal,1) > F,
  % repeat feature predictivity
  F = 17;
  l = -ones(size(s));
  l(label(s)==LABELS.intergenic) = +1;
  x = signal(F,s);
  [TP FP TN FN] = eval_separation(x, l, 10);
  OPT.line_style = {'.'};
  h = figure
  plot_ROC(TP, FP, TN, FN, OPT);
  title('repeat feature predictivity wrt. intergenic regions, ROC');
  grid on;
  % save figure as eps
  saveas(h,'repeat_predictivity.eps','eps');
end

%keyboard

% eof
