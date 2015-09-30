function show_gene_browser(CFG, test_chunks, anno, pred, cuffl, VIS_XCORR);


% normalize expression levels for mtim
if isfield(pred, 'expr'), 
  pred = norm_expr_levels(pred,10);
else
  for g=1:length(pred),
    pred(g).epxr = 10;
  end
end

% normalize expression levels for cufflinks if necessary (and existent)
if ~isempty(cuffl),
  if isfield(cuffl,'expr'),
    cuffl = norm_expr_levels(cuffl,10);
  else
    for g=1:length(cuffl),
      cuffl(g).epxr = 10;
    end
  end
end

% merge adjacent test chunks to save time during feature gathering
MAX_LEN = 20000;
c = size(test_chunks,1);
while c > 1,
  if test_chunks(c,3)-test_chunks(c,2)+1 > MAX_LEN,
    c = c-1;
    continue
  end

  % merge if test chunks are on the same chromosome and adjacent
  if (test_chunks(c,1)==test_chunks(c-1,1) && test_chunks(c,2)-1==test_chunks(c-1,3)),
    test_chunks(c-1,3) = test_chunks(c,3);
    % clear test chunk id
    test_chunks(c-1,4) = nan;
    test_chunks(c,:) = [];
  end
  c = c-1;
end
r = randperm(size(test_chunks,1));
test_chunks = test_chunks(r,:);

ag_int = [anno.start; anno.stop];
pg_int = [pred.start; pred.stop];
if ~isempty(cuffl),
  cg_int = [cuffl.start; cuffl.stop];
end

for c=1:size(test_chunks,1),
  chunk = test_chunks(c,:);
  signal = fill_chunks(chunk, CFG);
  offset = chunk(2);
  ag = anno([anno.chr_num]==chunk(1) ...
                & ag_int(2,:)>=chunk(2) & ag_int(1,:)<=chunk(3));
  pg = pred([pred.chr_num]==chunk(1) ...
                & pg_int(2,:)>=chunk(2) & pg_int(1,:)<=chunk(3));
  cg = [];
  if ~isempty(cuffl),
    cg = cuffl([cuffl.chr_num]==chunk(1) ...
                & cg_int(2,:)>=chunk(2) & cg_int(1,:)<=chunk(3));    
  end
  figure(1);
  plot_region(CFG, signal, offset, ag, pg, cg);

  % EVAL
%  fprintf('Evaluation:\n');
%  evaluate(ag,pg);
  

  % SHOW cross-correlation of features vs. signal
  if VIS_XCORR,
      % convert the annotated genes into a continuous signal
      fprintf('Converting annotated genes (%i) into signal...',length(ag));
      label = zeros(1,size(signal,2));
      for i=1:length(ag),
          for j=1:size(ag(i).exons,1),
              exons = ag(i).exons{j};
              for k=1:size(exons,1),
                  start = exons(k,1) - offset + 1;
                  stop = exons(k,2) - offset + 1;
                 
                  label(start:stop) = +1;
              end
          end
      end
      fprintf('Done!\n');

      % xcorr for all features
      figure(2);
      hold on;
      fprintf('Calculating cross-correlation coefficients..\n');
      
      MAX_LAG = length(label);
      %MAX_LAG = 250;
      
      xcorrs = [];
      names = get_feature_set_mTIM(); 
      for d=1:size(signal,1),
        %fprintf(' Dimension %2i..',d);
        % convert -/+inf into zeros
        signal(d,isinf(signal(d,:))) = 0.0;
        foo = xcorr(label,signal(d,:)); 
        xcorrs = [xcorrs; foo];

        % abs & normalization such that the signal is between 0..1
        foo = abs(foo)./(1.25*max(abs(foo)));
        %foo = foo./(2*foo(MAX_LAG));
        plot([1 length(foo)],[d d],':g');
        plot(1:length(foo),foo+d,'-k','linewidth',2); 

        % create a low-pass filtered signal for signal 10:13 (splice features)
        if (d>=10 && d<=13),
          %winsize = 100;
          % low-pass
          %foo = filter(ones(1,winsize)/winsize,1,signal(d,:));
          % high-pass 
          %foo = filter([ones(1,50) zeros(1,4) -ones(1,50)],1,signal(d,:));
          %foo = xcorr(label,foo);
          % high-pass on label
          foo = filter([-1 0 +1],1,label);
          foo = xcorr(foo,signal(d,:));
          
          foo = abs(foo)./(1.25*max(abs(foo)));
          plot(1:length(foo),foo+d,'-r','linewidth',1); 
        end

        names{d} = sprintf('%2i. %s',d,strrep(names{d},'_','-'));
      end
      fprintf(' Done!\n');

      % display zero lag
      plot([MAX_LAG MAX_LAG],[0 d+1],':b');

      hold off
      legend(names); 
  end

  keyboard
end


