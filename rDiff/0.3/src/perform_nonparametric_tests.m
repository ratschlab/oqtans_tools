function []=perform_nonparametric_tests(CFG,genes,variance_function_nonparametric_1, variance_function_nonparametric_2)


%Get the gene expression
fprintf('Loading gene expression\n')
if isempty(CFG.Counts_gene_expression)
    EXPR_TAB_FILENAME=[CFG.out_base 'Gene_expression.tab'];
else
    EXPR_TAB_FILENAME=CFG.Counts_gene_expression;
end

try
    Gene_expression=importdata(EXPR_TAB_FILENAME,'\t',1);
    Gene_expression=Gene_expression.data;
catch
    error(['Could not open: ' EXPR_TAB_FILENAME])
end
  

%Get the counts
fprintf('Loading nonparametric region counts\n')
if isempty(CFG.Counts_rDiff_nonparametric)
    IN_FILENAME=[CFG.out_base 'Nonparametric_region_counts.mat'];
    load(IN_FILENAME,'Counts_rDiff_nonparametric')
else
    IN_FILENAME=[CFG.out_base CFG.Counts_rDiff_nonparametric];
    load(IN_FILENAMEc,'Counts_rDiff_nonparametric') 
end



if CFG.use_rproc
   JB_NR=1;
   JOB_INFO = rproc_empty();
end

PAR.variance_function_nonparametric_1=variance_function_nonparametric_1;
PAR.variance_function_nonparametric_2=variance_function_nonparametric_2;
if 1==1
%%% Perform the test
% configuration
if not(CFG.use_rproc)
    fprintf('Performing nonparametric testing\n')
end
%define the splits of the genes for the jobs
idx=[(1:size(genes,2))',ceil((1:size(genes,2))*CFG.rproc_num_jobs/size(genes,2))']; 
% submit jobs to cluster

for i = 1:CFG.rproc_num_jobs
    PAR.genes = genes(idx(idx(:,2)==i,1)); 
    PAR.Counts_rDiff_nonparametric=Counts_rDiff_nonparametric(idx(idx(:,2)==i,1),:);
    PAR.Gene_expression=Gene_expression(idx(idx(:,2)==i,1),:);
    CFG.rproc_memreq = 5000;
    CFG.rproc_par.mem_req_resubmit = [5000 10000 32000];
    CFG.rproc_par.identifier = sprintf('Pnp.%i-',i);  
    CFG.outfile_prefix=[CFG.out_base_temp 'P_values_nonparametric_' num2str(i) '_of_' num2str(CFG.rproc_num_jobs)];
    PAR.CFG=CFG;
    if CFG.use_rproc
        fprintf(1, 'Submitting job %i to cluster\n',i);
        JOB_INFO(JB_NR) = rproc('get_nonparametric_tests_caller', PAR,CFG.rproc_memreq, CFG.rproc_par, CFG.rproc_time); 
        JB_NR=JB_NR+1;
    else
        get_nonparametric_tests_caller(PAR);
    end
end
if CFG.use_rproc
    [JOB_INFO num_crashed] = rproc_wait(JOB_INFO, 60, 1, -1);
end
end
% Get the test results

%%% Generate the output files
fprintf('Reading temporary results\n')
P_values_rDiff_nonparametric=ones(size(genes,2),1);
P_values_rDiff_mmd=ones(size(genes,2),1);


P_values_rDiff_nonparametric_error_flag=cell(size(genes,2),1);
P_values_rDiff_mmd_error_flag=cell(size(genes,2),1);
NAMES=cell(size(genes,2),1);
%Field containing the errors
ERRORS_NR=[];
idx=[(1:size(genes,2))',ceil((1:size(genes,2))*CFG.rproc_num_jobs/size(genes,2))']; 
% Iterate over the result files to load the data from the count files
for j = 1:CFG.rproc_num_jobs
    IN_FILENAME=[CFG.out_base_temp 'P_values_nonparametric_' num2str(j) '_of_' num2str(CFG.rproc_num_jobs)];
    IDX=idx(idx(:,2)==j,1);
    try
        load(IN_FILENAME)
        for k=1:length(IDX)
            if isempty(P_VALS{k,1}) %Gene was not tested for
                                          %some reason
                P_values_rDiff_nonparametric_error_flag{IDX(k)}='NOT_TESTED';
                P_values_rDiff_mmd_error_flag{IDX(k)}='NOT_TESTED';
            else
                if not(isempty(P_VALS{k,1}))
                    NAMES{IDX(k)}=P_VALS{k,1};
                end
                COUNTER=2;
                %Get the results from rDiff.mmd
                if CFG.perform_mmd
                    if not(isempty(P_VALS{k,COUNTER}))
                        P_values_rDiff_mmd(IDX(k))=P_VALS{k,COUNTER}{1};
                        if (isempty(P_VALS{k,COUNTER}{2}))
                            P_values_rDiff_mmd_error_flag{IDX(k)}='NOT_TESTED';
                        else
                            P_values_rDiff_mmd_error_flag{IDX(k)}='OK';
                        end
                    end
                    COUNTER=COUNTER+1;
                end
                
                %Get the results from rDiff.parametric
                if CFG.perform_nonparametric
                    if not(isempty(P_VALS{k,COUNTER}))
                        P_values_rDiff_nonparametric(IDX(k))=P_VALS{k,COUNTER}{1};
			if length(P_VALS{k,COUNTER})>1
			  if iscell(P_VALS{k,COUNTER}{2}) && length(P_VALS{k,COUNTER}{2}{3})>3
			    P_values_rDiff_nonparametric(IDX(k))=min(10*min(P_VALS{k,COUNTER}{2}{3})+max(P_VALS{k,COUNTER}{2}{3})*(10/(CFG.bootstraps+1)),1);
			  end
			end
                        if (isempty(P_VALS{k,COUNTER}{2}))
                            P_values_rDiff_nonparametric_error_flag{IDX(k)}='NOT_TESTED';
                        else
                            P_values_rDiff_nonparametric_error_flag{IDX(k)}='OK';
                        end
                    end;
                end
                
            end
        end
    catch
        for k=1:length(IDX)
            P_values_rDiff_nonparametric_error_flag{IDX(k)}='NOT_TESTED';
            P_values_rDiff_mmd_error_flag{IDX(k)}='NOT_TESTED';
        end
        warning(['There was a problem loading: ' IN_FILENAME ])
        ERRORS_NR=[ERRORS_NR;j];
    end
end
if not(isempty(ERRORS_NR))
    warning('There have been problems loading some of the parametric test result files');
end

fprintf('Writing output files\n')
%Generate P-value table for rDiff.nonparametric
if CFG.perform_nonparametric
%Open file handler
    P_TABLE_FNAME=[CFG.out_base 'P_values_rDiff_nonparametric.tab'];
    fid=fopen(P_TABLE_FNAME,'w');
    
    %print header
    fprintf(fid,'gene\tp-value\ttest-status\n');
    
    %print lines
    for j=1:size(genes,2)
        fprintf(fid,'%s',NAMES{j});
        fprintf(fid,'\t%f\t%s\n',P_values_rDiff_nonparametric(j),P_values_rDiff_nonparametric_error_flag{j});
    end
    %close file handler
    fclose(fid)
end



%Generate P-value table for rDiff.mmd
if CFG.perform_mmd
%Open file handler
    P_TABLE_FNAME=[CFG.out_base 'P_values_rDiff_mmd.tab'];
    fid=fopen(P_TABLE_FNAME,'w');
    
    %print header
    fprintf(fid,'gene\tp-value\ttest-status\n');
    
    %print lines
    for j=1:size(genes,2)
        fprintf(fid,'%s',NAMES{j});
        fprintf(fid,'\t%f\t%s\n',P_values_rDiff_mmd(j),P_values_rDiff_mmd_error_flag{j});
    end
    %close file handler
    fclose(fid)
end


