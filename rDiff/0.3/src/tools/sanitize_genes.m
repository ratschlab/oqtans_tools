function [genes]=sanitize_genes(genes,CFG)
%This function removes trnascript and genes which have a invalid
%structure and recomputes the splicegraph


%Mark genes with eronous exon definitions
RM_GENES_IDX=[]; %genes to keep
for i=1:size(genes,2)
  %remove transcripts which have a length smaller than
  %readlength
  RM_TR_IDX=[];
  START_MIN=inf;
  STOP_MAX=0;
  for j=1:size(genes(i).transcripts,2)
    if sum(genes(i).exons{j}(:,2)-genes(i).exons{j}(:,1))< CFG.sequenced_length
      RM_TR_IDX=[RM_TR_IDX,j];
    else
      START_MIN=min(START_MIN,genes(i).exons{j}(1,1));
      STOP_MAX=max(STOP_MAX,genes(i).exons{j}(end,2));
    end 
  end
  if ~isempty(RM_TR_IDX)
     genes(i).exons(RM_TR_IDX)=[];
    genes(i).transcripts(RM_TR_IDX)=[];
    genes(i).start=START_MIN;
    genes(i).stop=STOP_MAX;
  end
  if genes(i).start>START_MIN
      genes(i).start=START_MIN;
  end
  if genes(i).stop<STOP_MAX
      genes(i).stop=STOP_MAX;
  end
  
  
  %Check if exons are eronous
  CHECK=1;
  
  if size(genes(i).transcripts,2)==0
    CHECK=0;
  end
  
  for j=1:size(genes(i).transcripts,2)
    for k=1:(size(genes(i).exons{j},1)-1)
      if (genes(i).exons{j}(k,2)> genes(i).exons{j}(k+1,1))
	CHECK=0;
	break;
      end
    end
    
    if isempty(genes(i).exons{j})
	CHECK=0;
	break;
    end
    if CHECK==0
      break
    end
  end
  
  if genes(i).stop-genes(i).start<=CFG.sequenced_length
    CHECK=0;
  end
  if CHECK==0
    RM_GENES_IDX=[RM_GENES_IDX;i];
    genes(i).do_not_quant=1;
  else
      genes(i).do_not_quant=0;
  end
end
%genes(RM_GENES_IDX)=[];

% Create splicegraph 
for i=1:size(genes,2)
  gene=genes(i);
  ALL_EXO=[];
  
  for j=1:size(gene.exons,2)
    ALL_EXO=[ALL_EXO;gene.exons{j}];          
  end
  ALL_EXO=unique(ALL_EXO,'rows');
  GRAPH=zeros(size(ALL_EXO,1),size(ALL_EXO,1));
  for j=1:size(gene.exons,2)
    for k=1:(size(gene.exons{j},1)-1)
      [A,B]=intersect(ALL_EXO,gene.exons{j}(k,:),'rows');
      [A,C]=intersect(ALL_EXO,gene.exons{j}(k+1,:),'rows');
      GRAPH(B,C)=1;
    end
  end
  GRAPH=GRAPH+GRAPH';
  genes(i).splicegraph{1}=ALL_EXO';
  genes(i).splicegraph{2}=GRAPH;
end


  