function plot_region(CFG, obs_seq, offset, pred_genes, true_genes, pred_genes2)

OPT.colors = [0.3 0.3 0.8;
              0.5 0.5 1.0;
              0.4 0.4 1.0;
              0.3 0.8 0.3;
              0.5 1.0 0.5;
              0.0 0.8 0.0;
              0.0 0.8 0.0;
              0.0 0.0 0.8;
              0.0 0.0 0.8;
              0.5 0.8 0.5;
              0.5 0.8 0.5;
              0.5 0.5 0.8;
              0.5 0.5 0.8;
              0.5 0.1 0.1;
              0.2 0.2 0.2;
              0.5 0.5 0.8;
              0.5 0.5 0.8;
              0.5 0.1 0.1;
              0.2 0.2 0.2;
              0.5 0.5 0.5;
              1.0 0.2 0.2;
              0.2 1.0 0.2];

OPT.line_style = {'-', '-', '--', '-', '-', 'v', '^', 'v', '^',  'v', '^', '^', 'v', ...
                 '--', ':', '-', '-','--', ':', '-', '-','-'};

OPT.cutoffs = [0; -inf; -inf; 0; -inf; 0; 0; 0; 0; 0.85; 0.85; 0.85; 0.85; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf;-inf];

OPT.rescale = [1; 0.5; 0.2; 1; 0.5; 2; 2; 2; 2; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];

OPT.trunc = [0 inf;
             -inf inf;
             -inf inf;
             0 inf;
             -inf inf;
             0 inf;
             0 inf;
             0 inf;
             0 inf;
             0 1;
             0 1;
             0 1;
             0 1;
             0 1000;
             0 1000;
             0 1000;
             -inf inf;
             -inf inf;
             -inf inf;
             -inf inf;
             -inf inf;
             -inf inf];

OPT.baseline = [2; 0; 0; 2; 0; 0; 0; -4; -4; 0; 0; -2; -2; 2; 2; 2; 2; 2; 2; 2; 2; 2];


OPT.transform = {'log10', 'symlog10', '', 'log10', 'symlog10', 'log10', 'log10', ...
                 'log10','log10', '', '', '', '', 'log10', 'log10', 'log10', '', '', '', '', '', ''};

assert(size(OPT.colors,1) >= size(obs_seq,1));

assert(size(OPT.colors,1) >= size(obs_seq,1));
assert(length(OPT.line_style) >= size(obs_seq,1));

clf;
hold on;


for f=1:size(obs_seq,1),
  idx = find(~isnan(obs_seq(f,:)) & obs_seq(f,:) >= OPT.cutoffs(f));
  obs_seq(f,obs_seq(f,:)<OPT.trunc(f,1)) = OPT.trunc(f,1);
  obs_seq(f,obs_seq(f,:)>OPT.trunc(f,2)) = OPT.trunc(f,2);
  obs_seq(f,:) = OPT.rescale(f) * obs_seq(f,:);
  y = obs_seq(f,idx);
  if ~isempty(OPT.transform{f}),
    eval(sprintf('y = %s(obs_seq(f,idx));', OPT.transform{f}));
%    y(y<0 | isnan(y)) = 0;
  end
  y = y + OPT.baseline(f);

%  if (any(ismember([14,15,16],f))),
%    plot(idx+offset, y, OPT.line_style{f}, 'Color', OPT.colors(f,:), 'linewidth',2);
%  else
    plot(idx+offset, y, OPT.line_style{f}, 'Color', OPT.colors(f,:));  
%  end
   fprintf('f = %i\n',f);
end

  
comb_genes = {pred_genes};
gene_colors = {[0.1 0.8 0; 0.1 0 0.8]};
gene_y = -3.5;
gene_height = 0.2;
if exist('true_genes', 'var'),
  comb_genes{end+1} = true_genes;
  gene_colors{end+1} = [0.3 0.9 0; 0.3 0 0.9];
  gene_y(end+1) = gene_y(end)-8*gene_height;
end
if exist('pred_genes2', 'var'),
  comb_genes{end+1} = pred_genes2;
  gene_colors{end+1} = [0 0.6 0; 0 0 0.6];
  gene_y(end+1) = gene_y(end)-8*gene_height;
end

for h=1:length(comb_genes),
  genes = comb_genes{h};
  
  for g=1:length(genes),
    expr = 1;
    if isfield(genes, 'expr'),
      expr = mean(genes(g).expr) / 10;
      if expr>1, expr = 1; end
    end
    gh = gene_y(h);
    if length(genes(g).transcripts) > 8,
      warning('skipping transcript isoforms 9-%i', length(genes(g).transcripts));
    end
    for t=1:min(8, length(genes(g).transcripts)),
      gh = gh - gene_height;
      exons = genes(g).exons{t};
      num_exons = size(exons,1);
      introns = [exons(1:end-1,2)+1, exons(2:end,1)-1];
      num_introns = size(introns,1);
      % draw introns
      x = introns;
      switch genes(g).strand,
       case '+',
        y = (gh+0.25*gene_height)*ones(num_introns,2);
        for i=1:num_introns,
          plot(x(i,:), y(i,:), 'Color', gene_colors{h}(1,:)*expr, ...
               'LineWidth', 2);
        end
       case '-',
        y = (gh-0.25*gene_height)*ones(num_introns,2);
        for i=1:num_introns,
          plot(x(i,:), y(i,:), 'Color', gene_colors{h}(2,:)*expr, ...
               'LineWidth', 2);
        end
      end
      % draw exons
      x = [exons(:,2)'; exons(:,2)'; exons(:,1)'; exons(:,1)'];
      switch genes(g).strand,
       case '+',
        y = [(gh+0.5*gene_height)*ones(1,num_exons); ...
             gh*ones(1,num_exons); ...
             gh*ones(1,num_exons); ...
             (gh+0.5*gene_height)*ones(1,num_exons)];
        fill(x, y, gene_colors{h}(1,:)*expr);
       case '-',
        y = [gh*ones(1,num_exons); ...
             (gh-0.5*gene_height)*ones(1,num_exons); ...
             (gh-0.5*gene_height)*ones(1,num_exons); ...
             gh*ones(1,num_exons)];
        fill(x, y, gene_colors{h}(2,:)*expr);
      end
    end
  end
end
ylim([-9 7]);
% eof
