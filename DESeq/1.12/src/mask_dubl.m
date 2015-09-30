function [new_genes]=mask_dubl(genes,THRESH);

CHROMOSOMES={};
COUNTER=1;
for i=1:size(genes,2)
  CHROMOSOMES{COUNTER}=genes(i).chr;
  COUNTER=COUNTER+1;
end
CHROMOSOMES=unique(CHROMOSOMES);


INFO=zeros(size(genes,2),4);
for i=1:size(genes,2)
  CHR_VAL=0;
  for chr= 1:length(CHROMOSOMES)	
    if  strcmp(genes(i).chr,CHROMOSOMES(chr))
      CHR_VAL=chr;
    end
  end	
  INFO(i,:)=[i,genes(i).start,genes(i).stop, CHR_VAL];
end

COUNTER=1;
new_genes=genes;
for chr= 1:length(CHROMOSOMES)	
  GENES_ON_CHR=INFO(INFO(:,4)==chr,:);
  [TEMP,POS]=sort(GENES_ON_CHR(:,2));
  GENES_ON_CHR=GENES_ON_CHR(POS,:);
  STARTS=GENES_ON_CHR(:,2);
  STOPS=GENES_ON_CHR(:,3);
  for i=1:(size(GENES_ON_CHR,1))
    MIN_START=find(STOPS>=STARTS(i),1,'first');
    MAX_STOP=find(STARTS<=STOPS(i),1,'last');
    if MIN_START==i
      MIN_START=[];
    end
    if MAX_STOP==i
      MAX_STOP=[];
    end
    EXONS=[];
    if not (isempty(MIN_START))
      for CURR=MIN_START:(i-1)
	if(not(isempty(genes(GENES_ON_CHR(CURR,1)).transcripts)))
	  for tra=1:size(genes(GENES_ON_CHR(CURR,1)).transcripts,2)
	    if(not(isempty(genes(GENES_ON_CHR(CURR,1)).exons)))
	      EXONS=[EXONS;genes(GENES_ON_CHR(CURR,1)).exons{tra}];
	    else
	      EXONS=[EXONS;genes(GENES_ON_CHR(CURR,1)).start,genes(GENES_ON_CHR(CURR,1)).stop];
	    end	  
	  end
	else
	  EXONS=[EXONS;genes(GENES_ON_CHR(CURR,1)).start,genes(GENES_ON_CHR(CURR,1)).stop];
	end	  
      end
    end
    if not (isempty(MAX_STOP))
      for CURR=(i+1):MAX_STOP
	if(not(isempty(genes(GENES_ON_CHR(CURR,1)).transcripts)))
	  for tra=1:size(genes(GENES_ON_CHR(CURR,1)).transcripts,2)
	    if(not(isempty(genes(GENES_ON_CHR(CURR,1)).exons)))
	      EXONS=[EXONS;genes(GENES_ON_CHR(CURR,1)).exons{tra}];
	    else
	      EXONS=[EXONS;genes(GENES_ON_CHR(CURR,1)).start,genes(GENES_ON_CHR(CURR,1)).stop];	    
	    end
	  end
	else
	  EXONS=[EXONS;genes(GENES_ON_CHR(CURR,1)).start,genes(GENES_ON_CHR(CURR,1)).stop];
	end
	
      end
    end
    if  not (isempty([MAX_STOP,MIN_START]))
      EXONS=EXONS(EXONS(:,2)>=STARTS(i),:);
      EXONS=EXONS(EXONS(:,1)<=STOPS(i),:);
      new_genes(GENES_ON_CHR(i,1)).non_unique_regions=EXONS;    
    else
        new_genes(GENES_ON_CHR(i,1)).non_unique_regions=[];    
    end
  end
  COUNTER=COUNTER+1;
end
