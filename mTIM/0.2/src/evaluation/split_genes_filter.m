function genes = split_genes_filter(CFG,pred)
% A filter that checks if two adjacent genes have been 
% splitted although it is 'easy' to see that they belong together.
% 
% written by Nico Goernitz, Berlin Institute of Techology, 2013

warning('Split genes filter is a prototype!');
fprintf('Join adjacent genes that belong together.\n');

% this is a very small distance between adjacent genes 
% to be a bit more sure about joining 
MIN_DIST = 500;

% allow the intron span std to be 15% of the signal
MAX_STD_MUL = 0.15;

% [n x 2] vector of candidate indices 
gene_inds = [];


MEAN_EXPR = mean([pred.expr]);

% check each chromosome seperately
chrs = unique([pred.chr_num]);
for strand=['+','-'],
    for c=1:length(chrs),
        % get all genes for chromosome 'c'
        cinds = find([pred.chr_num]==chrs(c) & [pred.strand]==strand);

        % expr
        expr = [pred(cinds).expr];
        cinds = cinds(expr>(MEAN_EXPR*1.0));

   
        % sort the genes according to start position
        [startpos, inds] = sort([pred(cinds).start],'ascend');
        lens = [pred(cinds).stop] - [pred(cinds).start] + 1;
        lens = lens(inds);

        % mark all genes that have less than MINI_DIST distance
        dists = startpos(2:end) - (startpos(1:end-1)+lens(1:end-1));
        cand_inds = find(dists<MIN_DIST);

        if ~isempty(cand_inds),
            gene_inds = [gene_inds; cinds(inds(cand_inds))',cinds(inds(cand_inds+1))'];
        end
    end
end
fprintf('%i pairs of %i genes are very close.. \n',size(gene_inds,1),length(pred));


% there is no support for merging multiple genes yet
% hence, delete entries with multiple gene occurance
inds = unique(gene_inds);
for i=1:length(inds),
    ind = inds(i);
    foo = find(gene_inds(:,1)==ind | gene_inds(:,2)==ind);
    if length(foo)>1,
        gene_inds = gene_inds(setdiff([1:size(gene_inds,1)],foo(2:end)),:);
    end
end
fprintf('%i pairs remain after multiple occurance test. \n',size(gene_inds,1));


% fill output genes structure with all genes that
% remain untouched
inds = setdiff([1:length(pred)],unique(gene_inds));
genes = pred(inds);

% generate an intron signal from for all candidate gene pairs 
% and join them 

chunks = [];
startpos = [1];
for i=1:size(gene_inds,1),
    ind1 = gene_inds(i,1);
    ind2 = gene_inds(i,2);

    chr = pred(ind1).chr_num;
    start = pred(ind1).start;
    stop = pred(ind2).stop;
    len = stop-start+1;

    chunks = [chunks; chr,start,stop,-1,-1];
    startpos = [startpos, startpos(end)+len];
end
fprintf('Generate signal vector with length %i.\n',startpos(end)-1);

% feature 4 is intron span
signal = append_read_feats(chunks,zeros(0,startpos(end)-1),CFG);

% get the mean intron span of both genes and the split region
mean_span = mean(signal(4,signal(4,:)>0));
fprintf('Mean span is %f.\n',mean_span);

unsplitted = 0;
for i=1:size(gene_inds,1),
    ind1 = gene_inds(i,1);
    ind2 = gene_inds(i,2);
    start = pred(ind1).start;

    % get the mean intron span of the split region
    split_start = pred(ind1).stop - start + startpos(i);
    split_stop = pred(ind2).start - start + startpos(i);
    split_span = mean(signal(4,split_start:split_stop));
    split_span_std = std(signal(4,split_start:split_stop));

    max_std = split_span * 0.025;

    % if the 
    if (split_span>(mean_span*0.5) && split_span_std<max_std),
        % add the genes
        %  genes = [genes, pred(ind1), pred(ind2)];
        
        % combine the genes
        pred(ind1).stop = pred(ind2).stop;
        pred(ind1).expr = mean([pred(ind1).expr pred(ind2).expr]);
        %pred(ind1).exons{1} = [pred(ind1).exons{1}];
        %pred(ind1).exons{2} = [pred(ind2).exons{1}];
        pred(ind1).exons{1} = [pred(ind1).exons{1}; pred(ind2).exons{1}];
        
        %pred(ind1).transcripts{1} = pred(ind1).transcripts{1};
        %pred(ind1).transcripts{2} = pred(ind2).transcripts{1};
        pred(ind1).transcripts{1} = sprintf('mTIM:combined_%i_%i',ind1,ind2);
        
        genes = [genes, pred(ind1)];
        unsplitted = unsplitted + 1;

        fprintf('Unsplit gene with intron_span=%f  and std=%f  (%i).\n',split_span, split_span_std, unsplitted);
    else
        % leave them separate
        genes = [genes, pred(ind1), pred(ind2)];
    end


end

fprintf('%i genes have been joint.\n',unsplitted);
fprintf('Returning %i genes.\n',length(genes));



