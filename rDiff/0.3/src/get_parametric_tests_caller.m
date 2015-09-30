function []=get_parametric_tests_caller(PAR)


CFG = PAR.CFG;
genes = PAR.genes;
OUT_STR='';

% add paths
addpath(CFG.paths);
%load local variables
data_dir=CFG.data_dir;
OUT_STR=[];

variance_function_parametric_1=PAR.variance_function_parametric_1;
variance_function_parametric_2=PAR.variance_function_parametric_2;

Counts_rDiff_parametric=PAR.Counts_rDiff_parametric;
Gene_expression=PAR.Gene_expression;

%clear variabe PAR
clear PAR;

NUMBER_OF_TESTS_PER_GENE=(CFG.perform_parametric+CFG.perform_poisson);
if NUMBER_OF_TESTS_PER_GENE==0
    return
end
P_VALS=cell(size(genes,2),NUMBER_OF_TESTS_PER_GENE+12);


%iterate over genes
for i=1:size(genes,2)
  %TEMP_COUNT contains the counts for the current gene
  TEMP_COUNT=cell(1,3);
  gene = genes(i);
  
  
  OLD_OUT_STR=OUT_STR;
  OUT_STR=['Current gene: ' gene.name ' '];
  %print progress
  if CFG.use_rproc
      fprintf([OUT_STR '\n'])
  else
      % Erase old progress
      fprintf(repmat('\b',1,length(OLD_OUT_STR)));
      fprintf([OUT_STR])
  end
  
  %set default return values
  P_VALS{i,1}=gene.name;
  
  %check that the gene has exons defined
  if isempty(gene.exons)
    P_VALS{i,4}='Exons field empty in gene structure';
    continue;
  end

  %check that the gene is longer than the Reads. Otherwise the
  %definition of regions does not makes sense
  if gene.stop-gene.start<CFG.sequenced_length+3
    P_VALS{i,4}='Gene to short';
    continue;
  end
  %perform 
  COUNTER=2;
  if CFG.perform_parametric
      [PV, INFO]= rDiff_parametric(CFG,gene,Counts_rDiff_parametric(i,:),Gene_expression(i,:),variance_function_parametric_1, variance_function_parametric_2);
      P_VALS{i,COUNTER}={PV, INFO};
      COUNTER=COUNTER+1;
  end
      
  if CFG.perform_poisson
      [PV, INFO]= rDiff_poisson(CFG,gene,Counts_rDiff_parametric(i,:),Gene_expression(i,:));
      P_VALS{i,COUNTER}={PV, INFO};
  end
  
end

fprintf('\n')
%Save the p-values
OUT_FILENAME=[CFG.outfile_prefix '.mat'];
save(OUT_FILENAME,'P_VALS')

