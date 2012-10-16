function galaxy_rdiff_web(anno_dir, track1, track2, out_fname, test_meth)
% galaxy_rdiff_web(anno_dir, track1, track2, out_fname, test_meth)
%
% -- input --
% anno_dir: directory of genes
% track1: name of BAM file 1
% track2: name of BAM file 2
% out_fname: name of result file with p-values
% test_meth: test method ('poisson' or 'mmd')

%addpath('~/svn/tools/ngs/');  
%addpath('~/svn/projects/RNASeq_galaxy/rdiff.web/');  
samtools_dir = '/mnt/galaxyTools/tools/oqtans/sw_require/samtools/';

if ~exist(sprintf('%s.bai', track1), 'file')
  command = sprintf('%s./samtools index %s', samtools_dir, track1);
  [s m] = unix(command);
  if ~exist(sprintf('%s.bai', track1), 'file')
    sprintf('\nbai file for %s could not be created\n', track1); 
  end
end

if ~exist(sprintf('%s.bai', track2), 'file')
  command = sprintf('%s./samtools index %s', samtools_dir, track2);
  [s m] = unix(command);
  if ~exist(sprintf('%s.bai', track2), 'file')
    sprintf('\nbai file for %s could not be created\n', track2); 
  end
end

[ret, timedate] = unix('date');
timedate(timedate==sprintf('\n')) = [];
fprintf(1, '\n*** rDiff started with %s test %s *** \n\n', test_meth, timedate);

%%%% load genes %%%%%
load(sprintf('%s/genes.mat', anno_dir), 'genes');

%%%% rDiff %%%%%
fprintf(1, 'testing %i genes for differential expression...\n', length(genes));
p_values = eval(sprintf('%s_test(genes, track1, track2)', test_meth)); 
fprintf(1, 'done.\n');

%%%%% write to txt file %%%%%
write_rdiff_result(genes, p_values, out_fname, test_meth);

[ret, timedate] = unix('date');
timedate(timedate==sprintf('\n')) = [];
fprintf(1, '\n*** rDiff finished with %s test %s *** \n\n', test_meth, timedate);
