function []=get_read_counts(CFG,genes)

if CFG.use_rproc
   JB_NR=1;
   JOB_INFO = rproc_empty();
end

%%% Get the read counts
if CFG.estimate_gene_expression==1
    for RUN=1:size(CFG.BAM_FILES,2)
        % configuration
        CFG.curr_bamfile = CFG.BAM_FILES{RUN};
        if not(CFG.use_rproc)
            fprintf('Getting gene expression for: %s\n', CFG.curr_bamfile);
        end
        tic
        %define the splits of the genes for the jobs
        idx=[(1:size(genes,2))',ceil((1:size(genes,2))*CFG.rproc_num_jobs/size(genes,2))']; 
        % submit jobs to cluster
        for i = 1:CFG.rproc_num_jobs
            PAR.genes = genes(idx(idx(:,2)==i,1)); 
            CFG.rproc_memreq = 5000;
            CFG.rproc_par.mem_req_resubmit = [10000 20000 32000];
            CFG.rproc_par.identifier = sprintf('Exp.%i-',i);  
            CFG.outfile_prefix=[CFG.out_base_temp CFG.NAMES{RUN} '_' num2str(i) '_of_' num2str(CFG.rproc_num_jobs) '.mat'];
            PAR.CFG=CFG;
            if CFG.use_rproc
                fprintf(1, 'Submitting job %i to cluster\n',i);
                JOB_INFO(JB_NR) = rproc('get_reads_caller', PAR,CFG.rproc_memreq, CFG.rproc_par, CFG.rproc_time); 
                JB_NR=JB_NR+1;
            else
                get_reads_caller(PAR);
            end
        end
        toc   
    end
    if CFG.use_rproc
        [JOB_INFO num_crashed] = rproc_wait(JOB_INFO, 60, 1, -1);
    end
end

%%% Generate the output files
%load the count files
%load the count files
READS_PER_GENE=zeros(size(genes,2),size(CFG.BAM_FILES,2));
GENE_EXPR=zeros(size(genes,2),size(CFG.BAM_FILES,2));

Counts_rDiff_parametric=cell(size(genes,2),size(CFG.BAM_FILES,2));
Counts_rDiff_nonparametric=cell(size(genes,2),size(CFG.BAM_FILES,2));


%Field containing the errors
ERRORS_NR=[];
idx=[(1:size(genes,2))',ceil((1:size(genes,2))*CFG.rproc_num_jobs/size(genes,2))']; 
% Iterate over the result files to load the data from the count files
for RUN=1:size(CFG.BAM_FILES,2)
    for j = 1:CFG.rproc_num_jobs
        IN_FILENAME=[CFG.out_base_temp CFG.NAMES{RUN} '_' num2str(j) '_of_' num2str(CFG.rproc_num_jobs) '.mat'];
	IDX=idx(idx(:,2)==j,1);
        try
            load(IN_FILENAME)
            for k=1:length(IDX)
	        %Check wether COUNTS is empty
		if isempty(COUNTS{k})
		    continue
		end
                %Get the number of reads mapping to a gene
                if not(isempty(COUNTS{k}{2}))
                    READS_PER_GENE(IDX(k),RUN)=COUNTS{k}{2};
                end
                
                %get the number of nonalternative reads if
                %possible. Otherwise use the number of reads
                %mapping to a gene as gene expression
                if not(isempty(COUNTS{k}{3}))
                    GENE_EXPE(IDX(k),RUN)=COUNTS{k}{3};
                else
                    GENE_EXPE(IDX(k),RUN)=READS_PER_GENE(IDX(k),RUN);
                end
                
                %get the Counts for rDiff.parametric
                if not(isempty(COUNTS{k}{6}))
                    Counts_rDiff_parametric{IDX(k),RUN}=COUNTS{k}{6};
                end
                
                %get the counts for rDiff.nonparametric
                if not(isempty(COUNTS{k}{4}))
                    Counts_rDiff_nonparametric{IDX(k),RUN}=COUNTS{k}{4};
                end
            end
        catch
            warning(['There was a problem loading: ' IN_FILENAME ])
            % If something went wrong
            for k=1:length(IDX)
                READS_PER_GENE(IDX(k),RUN)=0;
                GENE_EXPR(IDX(k),RUN)=0;
                Counts_rDiff_parametric{IDX(k),RUN}={};
                Counts_rDiff_nonparametric{IDX(k),RUN}={};
            end
            ERRORS_NR=[ERRORS_NR; [RUN,i]];
        end
    end
end   
if not(isempty(ERRORS_NR))
    warning('There have been problems loading some of the raw count files');
end

%If less than 10 reads use abulute number of reads
GENE_EXPR(GENE_EXPR<10)=READS_PER_GENE(GENE_EXPR<10);

%Generate gene expression tables

%Open file handler for the gene expression table
EXPR_TAB_FILENAME=[CFG.out_base 'Gene_expression.tab'];
fid=fopen(EXPR_TAB_FILENAME,'w');

%print header
fprintf(fid,'gene');
for i=1:length(CFG.NAMES)
    fprintf(fid,'\t%s',CFG.NAMES{i});
end
fprintf(fid,'\n');

for j=1:size(genes,2)
   fprintf(fid,'%s',genes(j).name);
   for i=1:length(CFG.NAMES)
       fprintf(fid,'\t%i',GENE_EXPR(j,i));
   end
   fprintf(fid,'\n');
end

fclose(fid)

%Determine interpreter
if size(ver('Octave'),1)
    INTERPR = 1;
else
    INTERPR = 0;
end


%Save alternative region count file for rDiff.parametric
OUT_FILENAME=[CFG.out_base 'Alternative_region_counts.mat'];
if INTERPR
    save('-mat7-binary',OUT_FILENAME,'Counts_rDiff_parametric')
else
    save(OUT_FILENAME,'Counts_rDiff_parametric','-v7.3')
end

%Save alternative region count file for rDiff.nonparametric
OUT_FILENAME=[CFG.out_base 'Nonparametric_region_counts.mat'];
if INTERPR
    save('-mat7-binary',OUT_FILENAME,'Counts_rDiff_nonparametric')
else
    save(OUT_FILENAME,'Counts_rDiff_nonparametric','-v7.3')
end

return
