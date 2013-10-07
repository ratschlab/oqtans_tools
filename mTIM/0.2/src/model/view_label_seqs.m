function fh = view_label_seqs(fh, obs_seq, true_label_seq, pred_label_seq, second_label_seq)

% fh = view_label_seqs(fh, obs_seq, true_label_seq, [pred_label_seq], [second_label_seq])
%
% Visualizes a label sequence together with observation sequence and
% predicted label sequences if specified.
%
% fh -- the figure handle to plot into
% obs_seq -- the feature matrix (sequence of observations) of size m x n
%   where m is equal to the number of features and n the lengthe of the
%   training sequence
% true_label_seq -- the true label sequence to be learned
% pred_label_seq -- optional parameter: a predicted label sequence of
%   length n
% second_label_seq -- optional paramter: another label sequence
%   (e.g. corresponding to the max-margin violator) of length n
% returns the figure handle into which the label sequences were plotted
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008-2009

LABELS = get_label_set_mTIM();
  
if nargin<3, error('at least 3 arguments expected'); end

figure(fh)
clf
hold on
c = colormap;
c = c(round(linspace(1,64,size(obs_seq,1))),:);
glyphs = '....^^vv++++...';
for i=1:2:3%12%size(obs_seq,1),
  idx = find(~isnan(obs_seq(i,:)) & obs_seq(i,:) > -Inf & obs_seq(i,:) < Inf);
  x = obs_seq(i,idx);
  x = x - min(obs_seq(i,idx));
  x = x / max(obs_seq(i,idx));
  plot(idx, x, glyphs(i), 'Color', c(i,:));
end
for i=2:2:4%12%size(obs_seq,1),
  idx = find(~isnan(obs_seq(i,:)) & obs_seq(i,:) > -Inf & obs_seq(i,:) < Inf);
  x = obs_seq(i,idx);
  x = x - min(x);
  x = x / max(x);
  plot(idx, x+1.5, glyphs(i), 'Color', c(i,:));
end
for i=13:14%12%size(obs_seq,1),
  idx = find(~isnan(obs_seq(i,:)) & obs_seq(i,:) > -Inf & obs_seq(i,:) < Inf);
  x = obs_seq(i,idx);
  x = x - min(x);
  x = x / max(x);
  plot(idx, x+3, glyphs(i), 'Color', c(i,:));
end
truth = nan(size(true_label_seq));
ew_idx = find(true_label_seq==LABELS.exon_W);
truth(ew_idx) = 1;
ec_idx = find(true_label_seq==LABELS.exon_C);
truth(ec_idx) = -1;
iw_idx = find(true_label_seq==LABELS.intron_W);
truth(iw_idx) = 0.3;
ic_idx = find(true_label_seq==LABELS.intron_C);
truth(ic_idx) = -0.3;
ig_idx = find(true_label_seq==LABELS.intergenic);
truth(ig_idx) = 0;

idx = unique([ew_idx, ec_idx, iw_idx, ic_idx, ig_idx]);
%plot(idx, truth(idx)+0.1, '+', 'Color', [0 0.5 0]);
idx = find(true_label_seq==LABELS.ambiguous);
%plot(idx, zeros(1,length(idx))+0.1, 'd', 'Color', [0.7 0.7 0.7]);

color = [0.5 0 0];
view_segmentation(true_label_seq, LABELS, -0.5, color)

if exist('pred_label_seq', 'var'),
  pred = nan(size(pred_label_seq));
  ew_idx = find(pred_label_seq==LABELS.exon_W);
  pred(ew_idx) = 1;
  ec_idx = find(pred_label_seq==LABELS.exon_C);
  pred(ec_idx) = -1;
  iw_idx = find(pred_label_seq==LABELS.intron_W);
  pred(iw_idx) = 0.3;
  ic_idx = find(pred_label_seq==LABELS.intron_C);
  pred(ic_idx) = -0.3;
  ig_idx = find(pred_label_seq==LABELS.intergenic);
  pred(ig_idx) = 0;
  
  idx = unique([ew_idx, ec_idx, iw_idx, ic_idx, ig_idx]);
%  plot(idx, pred(idx)+0.2, '+', 'Color', [0.8 0 0]);

  color = [0 0 0.5];
  view_segmentation(pred_label_seq, LABELS, -1.0, color);
end

if exist('second_label_seq', 'var'),
  pred = nan(size(second_label_seq));
  ew_idx = find(second_label_seq==LABELS.exon_W);
  pred(ew_idx) = 1;
  ec_idx = find(second_label_seq==LABELS.exon_C);
  pred(ec_idx) = -1;
  iw_idx = find(second_label_seq==LABELS.intron_W);
  pred(iw_idx) = 0.3;
  ic_idx = find(second_label_seq==LABELS.intron_C);
  pred(ic_idx) = -0.3;
  ig_idx = find(second_label_seq==LABELS.intergenic);
  pred(ig_idx) = 0;
  
  idx = unique([ew_idx, ec_idx, iw_idx, ic_idx, ig_idx]);
%  plot(idx, pred(idx)+0.3, '+', 'Color', [0.2 0 0.8]);
end

%axis([0, length(truth)+1 -10 10]);

% eof
