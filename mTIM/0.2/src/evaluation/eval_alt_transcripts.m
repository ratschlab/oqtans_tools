%pred = load(['/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/' ...
%           'thaliana/lsl/mTim/output_cleave_12/genome_wide_predictions/genes_rquant_filtered.mat'])

pred = load(['/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/' ...
             'elegans/lsl/mTim_new/output_cleave_12/genome_wide_predictions/genes_rquant_filtered.mat'])

cnt_const = 0;
cnt_alt = 0;
num_alt_trans = [];
for g=1:length(pred.genes),
  if length(pred.genes(g).transcripts) > 1,
    cnt_alt = cnt_alt + 1;
    num_alt_trans(end+1) = length(pred.genes(g).transcripts);
  else
    cnt_const = cnt_const + 1;
  end
end