function genes = half_open_to_closed(genes)
% genes = half_open_to_closed(genes)

for j = 1:length(genes)
	for k = 1:length(genes(j).exons)
		if genes(j).strand=='+'
			genes(j).exons{k}(:,2) = genes(j).exons{k}(:,2)-1;
			if isfield(genes, 'cds_exons') && iscell(genes(j).cds_exons) && ~isempty(genes(j).cds_exons{k})
				genes(j).cds_exons{k}(:,2) = genes(j).cds_exons{k}(:,2)-1;
			end
			if isfield(genes, 'utr_3prime') && iscell(genes(j).utr_3prime) && length(genes(j).utr_3prime)>=k && ~isempty(genes(j).utr_3prime{k})
				genes(j).utr_3prime{k}(:,2) = genes(j).utr_3prime{k}(:,2)-1;
			end
			if isfield(genes, 'utr3_exons') && iscell(genes(j).utr3_exons) && length(genes(j).utr3_exons)>=k && ~isempty(genes(j).utr3_exons{k})
				genes(j).utr3_exons{k}(:,2) = genes(j).utr3_exons{k}(:,2)-1;
			end
			if isfield(genes, 'utr_5prime') && iscell(genes(j).utr_5prime)&&  length(genes(j).utr_5prime)>=k && ~isempty(genes(j).utr_5prime{k})
				genes(j).utr_5prime{k}(:,2) = genes(j).utr_5prime{k}(:,2)-1;
			end
			if isfield(genes, 'utr5_exons') && iscell(genes(j).utr5_exons)&&  length(genes(j).utr5_exons)>=k && ~isempty(genes(j).utr5_exons{k})
				genes(j).utr5_exons{k}(:,2) = genes(j).utr5_exons{k}(:,2)-1;
			end
		else
			genes(j).exons{k}(:,1) = genes(j).exons{k}(:,1)+1;
			if isfield(genes, 'cds_exons') && iscell(genes(j).cds_exons) && ~isempty(genes(j).cds_exons{k})
				genes(j).cds_exons{k}(:,1) = genes(j).cds_exons{k}(:,1)+1;
			end
			if isfield(genes, 'utr_3prime') && iscell(genes(j).utr_3prime) && length(genes(j).utr_3prime)>=k && ~isempty(genes(j).utr_3prime{k})
				genes(j).utr_3prime{k}(:,1) = genes(j).utr_3prime{k}(:,1)+1;
			end
			if isfield(genes, 'utr3_exons') && iscell(genes(j).utr3_exons) && length(genes(j).utr3_exons)>=k && ~isempty(genes(j).utr3_exons{k})
				genes(j).utr3_exons{k}(:,1) = genes(j).utr3_exons{k}(:,1)+1;
			end
			if isfield(genes, 'utr_5prime') && iscell(genes(j).utr_5prime) && length(genes(j).utr_5prime)>=k && ~isempty(genes(j).utr_5prime{k})
				genes(j).utr_5prime{k}(:,1) = genes(j).utr_5prime{k}(:,1)+1;
			end
			if isfield(genes, 'utr5_exons') && iscell(genes(j).utr5_exons) && length(genes(j).utr5_exons)>=k && ~isempty(genes(j).utr5_exons{k})
				genes(j).utr5_exons{k}(:,1) = genes(j).utr5_exons{k}(:,1)+1;
			end
		end
	end
end
