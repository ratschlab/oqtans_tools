function [COUNTS]=get_reads_caller(PAR)

CFG = PAR.CFG;
genes = PAR.genes;
clear PAR;

% add paths
addpath(CFG.paths);

data_dir=CFG.data_dir;
OUT_STR=[];

COUNTS=cell(size(genes,2),1);


% NON_PARAM_SAMPLE contains the read start density 
if CFG.perform_nonparametric
    NON_PARAM_SAMPLE=sparse([],[],[],10000,1,5000);
end


for i=1:size(genes,2)
  
    
      
  %TEMP_COUNT contins the counts for the current gene
  TEMP_COUNT=cell(1,7);
  gene = genes(i);
  
  %set default return values
  TEMP_COUNT{1}=gene.name;
  TEMP_COUNT{2}=[];
  TEMP_COUNT{3}=[];
  TEMP_COUNT{4}=[];
  TEMP_COUNT{5}=[];
  TEMP_COUNT{6}=[];
  
  %check that the gene has exons defined
  if isempty(gene.exons)
    STAT{i}=TEMP_COUNT;
    continue;
  end

  %check that the gene is longer than the Reads. Otherwise the
  %definition of regions does not makes sense
  if gene.stop-gene.start<CFG.sequenced_length+3
    STAT{i}=TEMP_COUNT;
    continue;
  end
  
  %Get the reads from gene
  [reads] = get_reads_for_gene(CFG,gene);

 
  %Get total number of reads
  NR_OF_READS=size(reads,1);
  TEMP_COUNT{2}=NR_OF_READS;
  
 
  NR_OF_TRANS=size(gene.transcripts,2);
  %check that the gene has more than one isoform
  if NR_OF_TRANS<=1
      STAT{i}=TEMP_COUNT;
      continue;
  end
 
  
  EXON_IDX=sum(gene.exonsequence,1)<NR_OF_TRANS;

  %Transform the reads in to counts per region
  [NEW_READS,UNEXPLAINED_READS,UNEXPLAINED_INDEX]= convert_reads_to_region_indicators(reads,gene);

  if length(unique(sum(NEW_READS,2)))>1
      warning(['Assignment of reads to regions is not unique for gene:' gene.name  '\n']); 
  end
  
  
  %Calulate gene expression
  %Get the non_unique_regions
  NON_ALT_EIRS=sum(gene.eirs_in_seq,1)==NR_OF_TRANS;
  TEMP_COUNT{3}=sum(sum(NEW_READS(:,NON_ALT_EIRS),2)>0); 
  
  
  %Get Counts for nonparametric variance function
  if CFG.perform_nonparametric
      %Get the read starts
      [TEMP,START]=max(reads,[],2); 
      read_starts=sparse((1:NR_OF_READS)',START,ones(NR_OF_READS,1),NR_OF_READS,size(reads,2),NR_OF_READS);
      
      %Get the index of the alternative reads
      ALT_READ_IX=zeros(size(reads,1),1);
      ALT_EIRS=and(sum(gene.eirs_in_seq,1)<NR_OF_TRANS,sum(gene.eirs_in_seq,1)>0);
      ALT_READ_IX(UNEXPLAINED_INDEX==0)=sum(NEW_READS(:,ALT_EIRS),2)>0;
      
      %Get the coverage of the read starts
      %COVERAGE=sum(read_starts(find(ALT_READ_IX>0),:),1);
      if CFG.only_gene_start
          COVERAGE=sum(reads,1);
      else
          COVERAGE=sum(reads(find(ALT_READ_IX>0),:),1);
      end
      if max(COVERAGE)>0
          TEMP_COUNT{4}=COVERAGE;
      else
          TEMP_COUNT{4}=[];
      end 
  end
  
  
  %Get counts for parametric settting
  if or(CFG.perform_parametric,CFG.perform_poisson)
      %Get the alternative reads 
      ALT_EIRS=and(sum(gene.eirs_in_seq,1)<NR_OF_TRANS,sum(gene.eirs_in_seq,1)>0);
      
      %Return the Counts in the alternative EIRS
      
      COUNTS_PER_EIR=sum(NEW_READS(:,ALT_EIRS),1);
      EXS_SEQ=gene.eirs_in_seq(:,ALT_EIRS);
      [NEWCOLS,IDX2,POS]=unique(EXS_SEQ','rows');
      NEWCOLS=NEWCOLS';
      EIR_COUNTS=zeros(1,length(IDX2));
      for j=1:max(POS)
          EIR_COUNTS(j)=sum(COUNTS_PER_EIR(POS==j));
      end
      TEMP_COUNT{6}=EIR_COUNTS;
  end
  
  clear NEW_READS
  clear reads;
  COUNTS{i}=TEMP_COUNT;
  
  OLD_OUT_STR=OUT_STR;
  OUT_STR=['Finished ' num2str(i) ' out of ' num2str(size(genes,2)) ' genes'];
  %print progress
  if CFG.use_rproc
      fprintf([OUT_STR '\n'])
  else
      % Erase old progress
      fprintf(repmat('\b',1,length(OLD_OUT_STR)));
      fprintf([OUT_STR])
  end
  
  
end
fprintf('\n')
%Save the counts
OUT_FILENAME=[CFG.outfile_prefix];
save(OUT_FILENAME,'COUNTS')
%Save the nonparametric histogram

