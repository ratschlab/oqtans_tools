function prepare_read_features(CFG)

% prepare_read_features(CFG)
%
% Generates positional features from short read alignment maps
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009


% TODO currently this is not used (porbably still buggy!), new code with
% corresponding functionality is inlined in make_chunks, should however at
% some point be in this function again
error('Don''t use this function!');

strands = '+-';

exon_read_map   = cell(1, CFG.num_chr);
intron_read_map = cell(2, CFG.num_chr);

% collect unspliced exons
for c=1:CFG.num_chr,
    %fprintf('    processing reads for %s \n', CFG.chr_names{c});
    pause(1);
    exon_read_map{c} = zeros(1, CFG.chr_lens(c), 'uint8');

    %Edited by Pramod
    region.start=1 ;
    region.stop=CFG.chr_lens(c) ;
    region.chr_num=c ;

    count=0; 
    for s=1:length(strands),
       region.strand = strands(s);
       fprintf('Processing reads for %s%s \n', CFG.chr_names{c}, region.strand);
       pause(1);
       [cover excluded_reads reads_ok] = get_coverage_per_read(CFG, region, 1);
       if(count==0),
          exon_read_map{c} = uint8(cover);
          count=1;
       else
          exon_read_map{c} = exon_read_map{c} + uint8(cover);
       end
       intron_read_map{s,c} = uint8(cover);
    end
    %End of edit
end 

% eof
