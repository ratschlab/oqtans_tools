function write_rdiff_result(genes, p_values, out_fname, test_meth) 
% function write_rdiff_result(genes, out_fname, source, mapped_reads, read_len) 
%
% -- input --
% genes: struct defining genes with start, stops, exons etc.
% p_values: p-values of test
% out_fname: name of result file with p-values
% test_meth: test method ('poisson' or 'mmd')

if exist(out_fname, 'file')
  fprintf('replacing file %s...\n', out_fname);
  [fd msg] = fopen(out_fname, 'w+');
  disp(msg);
else
  fprintf('creating file %s...\n', out_fname);
  [fd msg] = fopen(out_fname, 'w+');
  disp(msg);
end

fprintf(fd, '# gene name\ttesting method\tp-value\n');
for g = 1:length(genes),
  fprintf(fd, '%s\t%s\t%.4f\n', genes(g).name, test_meth, p_values(g));
end

fclose(fd);
