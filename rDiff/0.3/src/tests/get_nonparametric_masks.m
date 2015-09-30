function [MASKS]=get_nonparametric_masks(CFG,reads1,reads2)

%What fracitons should be choosen for the cutoff
cen_arr=0.1:0.1:1;

% Define the mask which should be used in order to mask high
% expresse genes
MASKS=zeros(length(cen_arr),size(reads1,2));

COUNTER=1;
for censor_frac= cen_arr  
    temp_reads1=reads1;
    temp_reads2=reads2;
    %cut to relvant position
    read_coverage=sum(reads1,1)+sum(reads2,1); 
    % get positions with a positive coverage
    nonzero_position=read_coverage>0;
    %Determine the cutoff values
    sorted_coverage=sort(read_coverage(nonzero_position));
    nr_of_nonzero_positions=sum(nonzero_position);
    relevant_positions=read_coverage<=sorted_coverage(ceil(nr_of_nonzero_positions*censor_frac));
    MASKS(COUNTER,relevant_positions)=1;
    COUNTER=COUNTER+1;
end
