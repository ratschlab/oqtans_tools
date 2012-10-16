function auc = plot_PRC(TP, FP, TN, FN, OPT)

% auc = plot_PRC(TP, FP, TN, FN, [OPT])
%
% Plots one or more precision recall curves and returns the AUC(s).
%
% OPT a struct with fields for the following options
%   colors -- optionally, colors for each graph can be specified
%   text_pos -- optionally, the position where the AUC is printed can be
%     specified as [pos_x, pos_y] with as many rows as TP, FP, TN, FN
%   ah -- optional axis handle to plot into
%   area_range -- restrict the calculated area under the curve to an
%     interval [area_range(1), area_range(2)]
%   line_style -- optional LineStyle property provided as cell array of
%     strings (with as many entries as TP, FP TN and FN have rows).
%
% returns the area under the precision-recall curve
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2007-2009

assert(size(TP,1) == size(FP,1));
assert(size(TP,1) == size(TN,1));
assert(size(TP,1) == size(FN,1));
num_curves = size(TP,1);

if exist('OPT', 'var') && isfield(OPT, 'text_pos'),
  assert(size(OPT.text_pos,1)==num_curves);
  assert(size(OPT.text_pos,2)==2);
end

if exist('OPT', 'var') && isfield(OPT, 'colors'),
  assert(size(OPT.colors,1)==num_curves);
  assert(size(OPT.colors,2)==3);
else
  c = color_grad();
  OPT.colors = c(round(linspace(1,64,num_curves)),:);
end

if exist('OPT', 'var') && isfield(OPT, 'line_style'),
  assert(length(OPT.line_style)==num_curves);
end

hold on
for i=1:num_curves
  if exist('OPT', 'var') && isfield(OPT, 'line_style'),
    line_style = OPT.line_style{i};
  else
    line_style = '-';
  end
  % Recall
  % add a tiny bit to denominator to avoid division by 0
  TPR = TP(i,:)./(TP(i,:)+FN(i,:) + 10^-6);
  if TP(i,1)==TP(i,1)+FN(i,1),
    TPR(1) = 1;
  end
  % Precision
  % add a tiny bit to denominator to avoid division by 0
  TDR = TP(i,:)./(TP(i,:)+FP(i,:) + 10^-6);
  if exist('OPT', 'var') && isfield(OPT, 'ah'),
    plot(OPT.ah, TPR, TDR, line_style, ...
         'Color', OPT.colors(i,:), 'LineWidth', 2);
  else
    plot(TPR, TDR, line_style, ...
         'Color', OPT.colors(i,:), 'LineWidth', 2);
  end
  
  % compute area under the curve
  TPR = TPR(end:-1:1);
  TDR = TDR(end:-1:1);
  assert(TPR(1)<0+10-5 & TPR(end)>1-10^-5);
  if exist('OPT', 'var') && isfield(OPT, 'area_range'),
    assert(TPR(1)<=OPT.area_range(1) & OPT.area_range(2)<=TPR(end));
    lb = find(TPR>OPT.area_range(1), 1, 'first');
    rb = find(TPR<OPT.area_range(2), 1, 'last');
    lalpha = (TPR(lb) - OPT.area_range(1)) / (TPR(lb) - TPR(lb-1));
    ralpha = (OPT.area_range(2) - TPR(rb)) / (TPR(rb+1) - TPR(rb));
    TDR_range = [TDR(lb)*(1-lalpha) + TDR(lb-1)*lalpha, ...
                 TDR(lb:rb), ...
                 TDR(rb)*(1-ralpha) + TDR(rb+1)*ralpha];
    TPR_range = [OPT.area_range(1), TPR(lb:rb), OPT.area_range(2)];
    auc = area_trapz(TPR_range, TDR_range);
    
    plot([OPT.area_range(1) OPT.area_range(1)], ...
         [0 TDR_range(1)], '--k');
    plot([OPT.area_range(2) OPT.area_range(2)], ...
         [0 TDR_range(end)], '--k');
  else
    auc = area_trapz(TPR, TDR);
  end
  if ~exist('OPT', 'var') || ~isfield(OPT, 'ah'),
    if ~exist('OPT', 'var') || ~isfield(OPT, 'text_pos'),
      th = text(0.1, 0.8-i*0.1, sprintf('AUC = %2.3f', auc));
    else
      th = text(OPT.text_pos(i,1), OPT.text_pos(i,2), ...
                sprintf('AUC = %2.3f', auc));
    end
    set(th, 'Color', OPT.colors(i,:));
  end
end

% determine the fraction of positive (labels)
% as a mesure of random guessing performance
for i=1:num_curves
  p = sum([TP(i,:); FN(i,:)]);
  p = unique(p);
  assert(length(p) == 1);
  n = sum([TN(i,:); FP(i,:)]);
  n = unique(n);
  assert(length(n) == 1);
  frac_pos(i) = p / (p + n);
end
frac_pos = unique(frac_pos);

if exist('OPT', 'var') && isfield(OPT, 'ah'),
  for f=1:length(frac_pos),
    plot(OPT.ah, [0 1], [frac_pos frac_pos], ':k');  
  end
  title(OPT.ah, 'Precision-recall curve');
  xlabel(OPT.ah, 'recall (sensitivity)');
  ylabel(OPT.ah, 'precision');
  axis(OPT.ah, [0 1 0 1]);
else
  for f=1:length(frac_pos),
    plot([0 1], [frac_pos frac_pos], ':k');  
  end
  title('Precision-recall curve');
  xlabel('recall (sensitivity)');
  ylabel('precision');
  axis([0 1 0 1]);
end

% eof