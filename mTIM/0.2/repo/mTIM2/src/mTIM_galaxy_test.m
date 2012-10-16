function mTIM_galaxy_test()

% test with drosophila
rgasp_dir = '/fml/ag-raetsch/nobackup/projects/rgasp.2/';
data_dir = '/fml/ag-raetsch/nobackup/projects/sequencing_runs/D_melanogaster/reads/';


% paths/info coming from galaxy
genome_dir = [rgasp_dir 'genomes/drosophila/drosophila.gio'];
anno_dir = [rgasp_dir 'annotations/drosophila/'];
splice_site_dir = [data_dir 'D_melanogaster_L3.6.bestpaired.filtered.bam.splice_pred/'];
acc_splice_dir = [splice_site_dir 'splice_acc.bspf/pred/'];
don_splice_dir = [splice_site_dir 'splice_don.bspf/pred/'];
read_map_file = [data_dir 'D_melanogaster_L3.6.bestpaired.filtered.bam'];
out_dir='/fml/ag-raetsch/home/ngoernitz/tmp/galaxy/testuser';

mTIM_galaxy_prepntrain(...
    genome_dir, acc_splice_dir, don_splice_dir, read_map_file, anno_dir, out_dir);




