function auc = plot_ROC(TP, FP, TN, FN, OPT)

% auc = plot_ROC(TP, FP, TN, FN, [OPT])
%
% Plots one or more ROC curves and returns the AUC(s).
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
% returns the area under the ROC curve
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2007-2008

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
for i=1:num_curves,
  if exist('OPT', 'var') && isfield(OPT, 'line_style'),
    line_style = OPT.line_style{i};
  else
    line_style = '-';
  end
  % 1-Specificity
  % add a tiny bit to denominator to avoid division by 0
  FPR = FP(i,:)./(FP(i,:)+TN(i,:) + 10^-6);
  if FP(i,1)==FP(i,1)+TN(i,1),
    FPR(1) = 1;
  end
  % Sensitivity (=Recall)
  % add a tiny bit to denominator to avoid division by 0
  TPR = TP(i,:)./(TP(i,:)+FN(i,:) + 10^-6);
  if exist('OPT', 'var') && isfield(OPT, 'ah'),
    plot(OPT.ah, FPR, TPR, line_style, 'Color', OPT.colors(i,:), 'LineWidth', 2);
  else
    plot(FPR, TPR, line_style, 'Color', OPT.colors(i,:), 'LineWidth', 2);
  end
  % compute area under the curve
  FPR = FPR(end:-1:1);
  TPR = TPR(end:-1:1);
  assert(FPR(1)<0+10-3 & FPR(end)>1-10^-3);
  if exist('OPT', 'var') && isfield(OPT, 'area_range'),
    assert(FPR(1)<=OPT.area_range(1) & OPT.area_range(2)<=FPR(end));
    lb = find(FPR>OPT.area_range(1), 1, 'first');
    rb = find(FPR<OPT.area_range(2), 1, 'last');
    lalpha = (FPR(lb) - OPT.area_range(1)) / (FPR(lb) - FPR(lb-1));
    ralpha = (OPT.area_range(2) - FPR(rb)) / (FPR(rb+1) - FPR(rb));
    TPR_range = [TPR(lb)*(1-lalpha) + TPR(lb-1)*lalpha, ...
                 TPR(lb:rb), ...
                 TPR(rb)*(1-ralpha) + TPR(rb+1)*ralpha];
    FPR_range = [OPT.area_range(1), FPR(lb:rb), OPT.area_range(2)];
    auc = area_trapz(FPR_range, TPR_range);
    
    plot([OPT.area_range(1) OPT.area_range(1)], ...
         [0 TPR_range(1)], '--k');
    plot([OPT.area_range(2) OPT.area_range(2)], ...
         [0 TPR_range(end)], '--k');
  else
    auc = area_trapz(FPR, TPR);
  end
  if ~exist('OPT', 'var') || ~isfield(OPT, 'ah'),
    if ~exist('OPT', 'var') || ~isfield(OPT, 'text_pos'),
      th = text(0.7, 0.8-i*0.1, sprintf('AUC = %2.3f', auc));
    else
      th = text(OPT.text_pos(i,1), OPT.text_pos(i,2), ...
                sprintf('AUC = %2.3f', auc));
    end
    set(th, 'Color', OPT.colors(i,:));
  end
end
if exist('OPT', 'var') && isfield(OPT, 'ah'),
  plot(OPT.ah, [0 1], [0 1], ':k');
  title(OPT.ah, 'ROC curve');
  xlabel(OPT.ah, '1 - specificity');
  ylabel(OPT.ah, 'sensitivity');
  axis(OPT.ah, [0 1 0 1]);
else
  plot([0 1], [0 1], ':k');
  title('ROC curve');
  xlabel('1 - specificity');
  ylabel('sensitivity');
  axis([0 1 0 1]);
end

% eof