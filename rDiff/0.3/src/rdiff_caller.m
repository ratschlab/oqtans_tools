function STAT = rDiff_caller(PAR)
% genes = opt_transcripts_caller(PAR)
%
% -- input --
% PAR contains
%     CFG: configuration struct
%     genes: struct defining genes with start, stops, exons etc.
%     profile_weights: weights of profile functions
%     intron_dists: distances to closest intron
%
% -- output --
% genes: struct with additional fields of eg. estimated transcript weights
addpath ../../tests
addpath ../../variance/final
addpath ../../experimental

samtools='/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/samtools/';
CFG = PAR.CFG;
SUFF_NAME=CFG.SUFF_NAME; 
%load(CFG.genes_path , 'genes'); 
%[genes]=mask_dubl(genes,0); 

load(CFG.genes_path , 'genes');
[genes]=mask_dubl(genes,0);

%[genes]=remove_transcripts(genes);
genes=genes(PAR.genes_idx);
BAM_FILES=CFG.BAM_FILES;
JB_NR=CFG.JB_NR;



%clear PAR;

%%%% paths
addpath(CFG.paths);
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/variance/nbin');
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/experiments/');
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/tests/');
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/read_utils/');
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/tests/sequence_tools/');
addpath(genpath('/fml/ag-raetsch/home/drewe/svn/tools/chronux/'))
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/tests/functions/');

out_base=CFG.out_base;
%sim_data_base=CFG.sim_data_base;
real_data_base=CFG.real_data_base;

expr12_bam=CFG.expr12_bam;
expr11_bam=CFG.expr11_bam;
expr22_bam=CFG.expr22_bam;
expr21_bam=CFG.expr21_bam;

CFG.SEQUENCED_LENGTH=80;
SEQUENCED_LENGTH=80;
%make map: gene->file

IDX1=CFG.IDX11;
IDX2=CFG.IDX22;

%BAM_FILES{IDX1}
%BAM_FILES{IDX2}


%genes=PAR.genes;
clear PAR

STAT=cell(size(genes));
do_save = 0;
out_base = '/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/out';
save_dir=fullfile(out_base,'cache');

for i=1:size(genes,2)
  %try
  RESULT=cell(1,7)
  gene = genes(i);
  
  RESULT{7}=JB_NR;
  RESULT{1}=gene.name
  load_only = false;
   
  if or(isempty(gene.exons),gene.stop-gene.start<=SEQUENCED_LENGTH)
    PV=1;
    RESULT{2}={PV,''};
    RESULT{3}={PV,''};
    RESULT{4}={PV,''};
    RESULT{5}=[Inf,Inf];
    RESULT{6}=[Inf,Inf];  
    RESULT{7}={'empty gene exons'};  
    RESULT{8}={PV,''};
    continue;
  end
  
  PV1=1;
  for l=1:size(gene.transcripts,2)
    EXONS=gene.exons{l}
    if sum(EXONS(:,2)-EXONS(:,1))<=SEQUENCED_LENGTH
          PV=1;
	  RESULT{2}={PV,''};
	  RESULT{3}={PV,''};
	  RESULT{4}={PV,''};
	  RESULT{5}=[Inf,Inf];
	  RESULT{6}=[Inf,Inf];  
	  RESULT{7}={'empty gene exons'};  
	  RESULT{8}={PV,''};
	  PV1=2;
	  
	  continue;
    end
  end
  if PV1==2
    continue;
  end
  gene.name 
    
   
  
  %adjust start and stop positions on the negative strand
  if strcmp(gene.strand,'-') 
    gene.start=gene.start+1;
    gene.stop=gene.stop+1;
    [reads11,reads12,intron11,intron12] = get_reads_dual2_intron(expr11_bam,expr12_bam,gene,samtools,false,10);
    [reads11,FLAG]=remove_reads_from_other_genes(reads11,gene);
    [reads12,FLAG]=remove_reads_from_other_genes(reads12,gene);
    
   
    [reads21,reads22,intron21,intron22] = get_reads_dual2_intron(expr21_bam,expr22_bam,gene,samtools,false,10);
    [reads21,FLAG]=remove_reads_from_other_genes(reads21,gene);
    [reads22,FLAG]=remove_reads_from_other_genes(reads22,gene);
    gene.start=gene.start-1;
    gene.stop=gene.stop-1;
  else
    [reads11,reads12,intron11,intron12] = get_reads_dual2_intron(expr11_bam,expr12_bam,gene,samtools,false,10);
    [reads11,FLAG]=remove_reads_from_other_genes(reads11,gene);
    [reads12,FLAG]=remove_reads_from_other_genes(reads12,gene);
    
    [reads21,reads22,intron21,intron22] = get_reads_dual2_intron(expr21_bam,expr22_bam,gene,samtools,false,10);
    [reads21,FLAG]=remove_reads_from_other_genes(reads21,gene);
    [reads22,FLAG]=remove_reads_from_other_genes(reads22,gene);
  end
  
  reads11=reads11(sum(reads11,2)>70,:);
  reads12=reads12(sum(reads12,2)>70,:);
  reads21=reads21(sum(reads21,2)>70,:);
  reads22=reads22(sum(reads22,2)>70,:);
  %% TEST with read cliping
  
  READS1={reads11,reads12};
  READS2={reads21,reads22};
  reads1=[reads11;reads12]; 
  reads2=[reads21;reads22]; 
	  
  TOTAL_SIZE=(size(reads12,1)+size(reads12,1)+(size(reads21,1)+size(reads22,1)));
  MIN_SIZE=min(size(reads12,1)+size(reads12,1),(size(reads21,1)+size(reads22,1)));
  %keyboard
  
  CLIP=0;
  if 1==0
  if (size(reads1,1)>0)
    for i=1:size(reads1,1)
      reads1(i,find(reads1(i,:),2,'first'))=0;
    end
    for i=1:size(reads1,1)
      reads1(i,find(reads1(i,:),2,'last'))=0;
    end
  end
  if (size(reads2,1)>0)
    for i=1:size(reads2,1)
      reads2(i,find(reads2(i,:),2,'first'))=0;
    end    
    for i=1:size(reads2,1)
      reads2(i,find(reads2(i,:),2,'last'))=0;
    end
  end
  CLIP=4;
  end
  
 % keybloard
  %% END TEST
  
  COUNTER=3;  
  FLAG=0;
  for k=1:size(gene.transcripts,2)
    if sum(gene.exons{k}(:,2)-gene.exons{k}(:,1))<75
      FLAG=1;
    end
  end

  if FLAG==1
    RESULT{2}={'gene shorter than 75 bp'};
    STAT{i}=RESULT;
    continue
  end
 
  if 1==0
      
      for V_COUNTER=2%:length(CFG.VARIANCES1)
          if (TOTAL_SIZE>0)
              if(MIN_SIZE>0)
                  PV=1;
                  INFO='';
                  [PV,INFO] =diff_mmd([reads11;reads12],[reads21;reads22],gene);
                  PV
                  RESULT{COUNTER}={PV,INFO};
                  COUNTER=COUNTER+1;
              else
                  PV=1;
                  RESULT{COUNTER}={PV,2};
                  COUNTER=COUNTER+1;
              end
          else
              PV=1;
              RESULT{COUNTER}={PV,1};
              COUNTER=COUNTER+1;
          end
      end 
  end
  % keyboard
  
  VAR_FACT=1;
  COV_POS=1;
  REWEIGHT=0;
  if 1==1
  for V_COUNTER=2%:length(CFG.VARIANCES1)
      if (TOTAL_SIZE>0)
          if(MIN_SIZE>0)
              cen_ARR=0.1:0.1:1;
              REWEIGHT_WEIGHTS=zeros(length( cen_ARR),size(reads1,2));
              cen_count=1;
              for censor_frac= cen_ARR
                  %if strcmp(gene.name,'AT1G01073')
                  %    keyboard
                  %end 
                  censor_frac
                  temp_reads1=reads1;
                  temp_reads2=reads2;
                  %cut to relevant positions
                  COVERAGE=sum(temp_reads1,1)+sum(temp_reads2,1); 
                  NONZERO=COVERAGE>0;
                  SORTED_COVERAGES=sort(COVERAGE(NONZERO));
                  NR_OF_NONZERO=sum(NONZERO);
                  CHOSEN_POSITIONS=COVERAGE<=SORTED_COVERAGES(ceil(NR_OF_NONZERO*censor_frac));
                  REWEIGHT_WEIGHTS(cen_count,CHOSEN_POSITIONS)=1;
                  cen_count=cen_count+1;
                  %temp_reads1=temp_reads1(:,CHOSEN_POSITIONS);
                  %temp_reads2=temp_reads2(:,CHOSEN_POSITIONS);
                  %remove reads which have no coverage;
                  %temp_reads1=temp_reads1(sum(temp_reads1,2)>0,:);
                  %temp_reads2=temp_reads2(sum(temp_reads2,2)>0,:);
                  %SIZE1=size(temp_reads1,1);
                  %SIZE2=size(temp_reads2,1);
                  %MIN_SIZE=min(SIZE1,SIZE2);
                  %TOTAL_SIZE=SIZE1+SIZE2;
              end
              if (TOTAL_SIZE>0)
                  if(MIN_SIZE>0)
                      PV=1;
                      INFO='';
                      REWEIGHT=1;
                      [PV,INFO] =diff_mmd_variance_subsample3({reads11,reads12},{reads21,reads22},gene,CFG.VARIANCES1{V_COUNTER},CFG.VARIANCES2{V_COUNTER},1,1,REWEIGHT, REWEIGHT_WEIGHTS)
                      RESULT{COUNTER}={PV,INFO};
                      COUNTER=COUNTER+1;
                  else
                      PV=1;
                      RESULT{COUNTER}={PV,2};
                      COUNTER=COUNTER+1;
                  end
              else
                  PV=1;
                      RESULT{COUNTER}={PV,1};
                      COUNTER=COUNTER+1;
              end
          end 
      end
  end
  end

  if 1==0  
    VAR_FACT=1;
    COV_POS=1;
    REWEIGHT=0;
    for V_COUNTER=2%:length(CFG.VARIANCES1)
      for cov_fact_ix=1:length(COV_POS)
	
	if (TOTAL_SIZE>0)
	  if(MIN_SIZE>0)
	    PV=1;
	    INFO='';
	    [PV,INFO] =diff_mmd_variance_subsample2({reads11,reads12},{reads21,reads22},gene,CFG.VARIANCES1{V_COUNTER},CFG.VARIANCES2{V_COUNTER},VAR_FACT(var_fact_ix),COV_POS,REWEIGHT);
	    RESULT{COUNTER}={PV,INFO};
	    COUNTER=COUNTER+1;
	  else
	    PV=1;
	    RESULT{COUNTER}={PV,2};
	    COUNTER=COUNTER+1;
	  end
	else
	  PV=1;
	  RESULT{COUNTER}={PV,1};
	  COUNTER=COUNTER+1;
	end
      end 
    end
      
    if 1==1
      VAR_FACT=1;
      COV_POS=1;
      REWEIGHT=1;
      for V_COUNTER=2%:length(CFG.VARIANCES1)
	for var_fact_ix=1:length(VAR_FACT)
	  
	  if (TOTAL_SIZE>0)
	    if(MIN_SIZE>0)
	      PV=1;
	      INFO='';
	      [PV,INFO] =diff_mmd_variance_subsample2({reads11,reads12},{reads21,reads22},gene,CFG.VARIANCES1{V_COUNTER},CFG.VARIANCES2{V_COUNTER},VAR_FACT(var_fact_ix),COV_POS,REWEIGHT)
	      RESULT{COUNTER}={PV,INFO};
	      COUNTER=COUNTER+1;
	    else
	      PV=1;
	      RESULT{COUNTER}={PV,2};
	      COUNTER=COUNTER+1;
	    end
	  else
	    PV=1;
	    RESULT{COUNTER}={PV,1};
	    COUNTER=COUNTER+1;
	  end
	end 
      end
    end
  end
  
  %keyboard
  TOTAL_SIZE=(size(reads12,1)+size(reads12,1)+(size(reads21,1)+size(reads22,1)));
  MIN_SIZE=min(size(reads12,1)+size(reads12,1),(size(reads21,1)+size(reads22,1)));
  if 1==0
    for V_COUNTER=1
      if (TOTAL_SIZE>0)
	if(MIN_SIZE>0)
	  PV=1;
	  INFO='';
	  if isempty([intron11,intron12])|isempty([intron21,intron22])
	    [PV,INFO] =diff_mmd_variance_NB_NB_simple([reads11;reads12],[reads21;reads22],gene,CFG.VARIANCES1{V_COUNTER},CFG.VARIANCES2{V_COUNTER});  
	  else	  
	    [PV,INFO] =diff_mmd_variance_splice([reads11;reads12],[reads21;reads22],0.5,[intron11,intron12],[intron21,intron22],gene,CFG.VARIANCES1{V_COUNTER},CFG.VARIANCES2{V_COUNTER});
	  end
	  RESULT{COUNTER}={PV,INFO};
	  COUNTER=COUNTER+1;
	else 
	  PV=1;
	  RESULT{COUNTER}={PV,2};
	  COUNTER=COUNTER+1;
	end
      else
	PV=1;
	RESULT{COUNTER}={PV,1};
	COUNTER=COUNTER+1;
      end
    end 
  end
  PV=1;
  %keyboard
  for V_COUNTER=length(CFG.VARIANCES1)
    if (TOTAL_SIZE>0)

      if(MIN_SIZE>0)   
	[P_VALUE, INFO]= diff_nbin7(READS1,READS2,gene,SEQUENCED_LENGTH,CFG.VARIANCES1{V_COUNTER},CFG.VARIANCES2{V_COUNTER});
  	RESULT{COUNTER}={P_VALUE,INFO};
	if not(isempty(INFO))
	  if iscell(INFO)
	    RESULT{COUNTER}=INFO{5};
	    COUNTER=COUNTER+1;
	  end
	end
      else
	PV=1;
	RESULT{COUNTER}={PV,2};
	COUNTER=COUNTER+1;
      end
    else
      PV=1;
      RESULT{COUNTER}={PV,1};
      COUNTER=COUNTER+1;
    end
  end
  

  [SPLICINGEVENTS,SEQUENCE,EXONSEQUENCE]=splicingsequence(gene);
  [UNIQUE_NEW_EXONS,GRAPHNODES,ORDER_OF_GRAPHNODE,EIRS_IN_SEQ]=transform_single_end_reads(SPLICINGEVENTS,SEQUENCE,EXONSEQUENCE,CFG.SEQUENCED_LENGTH-CLIP);
  [NEW_READS1,UNEXPLAINED_READS1,UNEXPLAINED_INDEX1]= convert_reads_to_region_indicators([reads11;reads12],UNIQUE_NEW_EXONS,GRAPHNODES,ORDER_OF_GRAPHNODE,EIRS_IN_SEQ,gene);
  [NEW_READS2,UNEXPLAINED_READS2,UNEXPLAINED_INDEX2]= convert_reads_to_region_indicators([reads21;reads22],UNIQUE_NEW_EXONS,GRAPHNODES,ORDER_OF_GRAPHNODE,EIRS_IN_SEQ,gene);
 % keyboard
  TOTAL_SIZE=(size(NEW_READS1,1)+(size(NEW_READS2,1)));
  MIN_SIZE=min(size(NEW_READS1,1),(size(NEW_READS2,1)));

  if (TOTAL_SIZE>0)
    if(MIN_SIZE>0)
      [PV,INFO] =  diff_poisson_bonf_3_unequal_segment(NEW_READS1,NEW_READS2,gene,CFG.SEQUENCED_LENGTH-CLIP);
      RESULT{COUNTER}={PV,INFO};
      COUNTER=COUNTER+1;
    else
      PV=1;
      RESULT{COUNTER}={PV,2};
      COUNTER=COUNTER+1;
    end
  else
    PV=1;
    RESULT{COUNTER}={PV,1};
    COUNTER=COUNTER+1;
  end
  PV=1
  if 1==0
  if (TOTAL_SIZE>0)
      if(MIN_SIZE>0)
      [PV,INFO] =  diff_poisson_bonf_4_unequal_segment(NEW_READS1,NEW_READS2,gene,CFG.SEQUENCED_LENGTH-CLIP);
      RESULT{COUNTER}={PV,INFO};
      COUNTER=COUNTER+1;
    else
      PV=1;
      RESULT{COUNTER}={PV,2};
      COUNTER=COUNTER+1;
    end
  else
    PV=1;
    RESULT{COUNTER}={PV,1};
    COUNTER=COUNTER+1;
  end
  %keyboard
  end
  STAT{i}=RESULT;
  
end;
%keyboard
OUT_FILENAME=['/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/out/analysis_artificial_variance_03_08_2012/' SUFF_NAME '_rep_mmd_07_07_2012_' int2str(JB_NR)   '.mat'];
save(OUT_FILENAME,'STAT')
