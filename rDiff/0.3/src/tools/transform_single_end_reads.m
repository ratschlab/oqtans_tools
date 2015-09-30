function [UNIQUE_NEW_EXONS, GRAPHNODES, ORDER_OF_GRAPHNODE, EIRS_IN_SEQ] = transform_single_end_reads(gene, SEQUENCED_LENGTH)
% This function calculates all regions onto which a read may fall. 
% SEQUENCED_LENGTH is the length of a read.

    SPLICINGEVENTS=gene.splicingevents;
    SEQUENCE=gene.sequence;
    EXONSEQUENCE=gene.exonsequence;
    
    NB_OF_TRANS = size(EXONSEQUENCE,1);
    
    NEWEXONS = 1:(length(SPLICINGEVENTS) - 1);
    UNIQUE_NEW_EXONS = NEWEXONS(SEQUENCE>0);
    
    NB_OF_EXONS = length(UNIQUE_NEW_EXONS);
    MAX_EXON_NB_TO_POS = cumsum((NEWEXONS .* SEQUENCE) > 0);
 
    % NB_EXONSEQUENCE is for each transcript the position in
    % SPLICINGEVENTS where one of its exonic region starts and 0 if it
    % is intronic
    
    NB_EXONSEQUENCE = EXONSEQUENCE .* repmat(NEWEXONS,NB_OF_TRANS,1);
    
    %%% This contains all exons
    GRAPHNODES = []; 
    
    %%% This is the array where a 1 is in column j if that EIR (Exons in
    %a region covered by a read) is contained in the transcript 
    EIRS_IN_SEQ = [];
        
    CURRENT_NODE = 0;% To catch errors with non-initalized variable		  
    
    EIR_TRANSCRIPTS = cell(1,NB_OF_TRANS);
    ORDER_OF_GRAPHNODE = cell(1,NB_OF_TRANS);
    LENGTHS_OF_GRAPHNODE = cell(1,NB_OF_TRANS);
    
    SPLICING_LENGTHS = SPLICINGEVENTS(2:end) - SPLICINGEVENTS(1:end-1);
    
    for i = 1:NB_OF_TRANS
        
        %%% remove zero positions, adjust splicing positions
        CURRENT_EXONS = NB_EXONSEQUENCE(i, EXONSEQUENCE(i,:) > 0);
        
        SPLICING_CORRECT = cumsum(SPLICING_LENGTHS .* (NB_EXONSEQUENCE(i,:) == 0));
        %SPLICING_CORRECT contains the position of SPLICINGSEQUENCE in
        %the transcript when all introns are spliced out
        SPLICING_CORRECT = [1, SPLICINGEVENTS(2:end) - SPLICING_CORRECT];
        
        %This ensures that the end of the transcript is also in CURRENT_SPLICINGEVENTS
        IDX = EXONSEQUENCE(i,:) == 1 ;
        IDX(find(EXONSEQUENCE(i,:) == 1, 1, 'last') + 1) = true;
        CURRENT_SPLICINGEVENTS = SPLICING_CORRECT(IDX);
        if length(CURRENT_SPLICINGEVENTS) == 0
            continue
        end
        LASTPOS = CURRENT_SPLICINGEVENTS(end);
        if LASTPOS <= SEQUENCED_LENGTH
            warning('CURRENT_SPLICINGEVENTS(end) > SEQUENCED_LENGTH') 
        end
        %assert(LASTPOS > SEQUENCED_LENGTH,'CURRENT_SPLICINGEVENTS(end) > SEQUENCED_LENGTH')
        % Calculate the positions when the EIRS can change
        
        % Determine positions which start SEQUENCED_LENGTH positions before a splicing event
        % defines a window of size SEQUENCED_LENGTH around the CURRENT_SPLICINGEVENTS
        READEVENTS_START = max([CURRENT_SPLICINGEVENTS(1:end - 1) - SEQUENCED_LENGTH + 1; ones(1,length(CURRENT_SPLICINGEVENTS)-1)],[],1);
        READEVENTS_END = min([CURRENT_SPLICINGEVENTS(2:end); repmat(LASTPOS - SEQUENCED_LENGTH,1,length(CURRENT_SPLICINGEVENTS(2:end)))],[],1);
        
        % Calculate EIRS 
        % CHANGE_POINTS are those points in a transcript where a EIR changes, namly the splicesites of that transcript plus and
        % minus the SEQUENCED_LENGTH - the above descibed window
        CHANGE_POINTS = unique([READEVENTS_START, READEVENTS_END]);
        
        for j = 1:(length(CHANGE_POINTS) - 1)
            
            POINTS_OF_INTEREST = ( READEVENTS_START(1,:) <= CHANGE_POINTS(j)) & (READEVENTS_END(1,:) > CHANGE_POINTS(j));
            
            % MAX_EXON_NB_TO_POS is mapping back to the unspliced coordinates
            CURRENT_EIRS = zeros(1,NB_OF_EXONS);
            CURRENT_EIRS( MAX_EXON_NB_TO_POS(CURRENT_EXONS(POINTS_OF_INTEREST))) = 1;
            
            %%% Already seen such exon composition in sliding
            %%% window?
            [TEMP, CURRENT_NODE] = intersect(GRAPHNODES, CURRENT_EIRS, 'rows');
            if isempty(TEMP)
		GRAPHNODES = [GRAPHNODES; CURRENT_EIRS];  %Add Key
		EIRS_IN_SEQ = [EIRS_IN_SEQ, zeros(NB_OF_TRANS,1)];
		CURRENT_NODE = size(GRAPHNODES,1);
            end
            
            EIRS_IN_SEQ(i,CURRENT_NODE) = 1;
            ORDER_OF_GRAPHNODE{i} = [ORDER_OF_GRAPHNODE{i}, CURRENT_NODE];
            LENGTHS_OF_GRAPHNODE{i} = [LENGTHS_OF_GRAPHNODE{i}, [CHANGE_POINTS(j); CHANGE_POINTS(j+1)]]; 
        end
    end
    