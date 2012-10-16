function [str,valid_splice,ok] = load_genomic(contig,strand,start,stop, genome_info,check_splice)
% function [str,valid_splice,ok]=load_genomic(contig,strand,start,stop, genome_info,check_splice)

if nargin<6
  check_splice=0;
end ;
ok=1 ;

if ~check_splice 
  valid_splice=[];
  % fprintf(1,'no splice site consensus check!')
end ;


if length(start)>1 | length(stop)>1,
  if nargout<3
    assert(length(start)==length(stop)) ;
    assert(all([start(2:end)-stop(1:end-1)]>=0));
  else
    if ~(length(start)==length(stop)) || ~(all([start(2:end)-stop(1:end-1)]>=0)),
      ok=0 ;
      str='' ;
      valid_splice = [] ;
      return ;
    end ;
  end ;
  valid_splice = zeros(length(start),2); 
  if strand == '+', 
    idx = 1:length(start) ; 
  else 
    idx = length(start):-1:1 ; 
  end ;
  str = '' ;
  for i=idx  %1:length(start)
    seq = load_genomic(contig, strand, start(i)-5, stop(i)+5, genome_info) ;
    str = [str seq(6:end-5)] ;

    % check internal splicing sites
    if check_splice
      % if (i~=1 && strand == '+') || (i~=length(start) && strand == '-')
      if (strand == '+') || (i~=length(start) && strand == '-')
        valid_splice(i,1) = isequal(seq(4:5),'ag') ;
        %if ~isequal(seq(4:5),'ag'), fprintf('acc: %s %c\n', seq(3:6), strand); end ;
        % assert(isequal(seq(4:5),'ag'))
      end ;
      %if i~=idx(end)
      % if (i~=1 && strand == '-') || (i~=length(start) && strand == '+')
      if (strand == '-') || (i~=length(start) && strand == '+')        
        valid_splice(i,2) = ismember(seq(end-5:end-4), {'gt','gc'}) ;
        %if ~isequal(seq(end-4:end-3),'gt')||isequal(seq(end-4:end-3),'gc'), fprintf('don: %s %c\n', seq(end-5:end-2), strand); end ;
        % assert(isequal(seq(end-4:end-3),'gt')||isequal(seq(end-4:end-3),'gc'))
      end ; 
    end ;
  end
  return ;
else
  valid_splice=[0 0];
end ;


if ischar(contig),
  %contig_idx = find(ismember(upper(genome_info.contig_names),upper(contig))) ;
  contig_idx = strmatch(upper(contig), upper(genome_info.contig_names), 'exact') ;
  assert(all(size(contig_idx)==1) & ~isempty(contig_idx)) ;
else
  contig_idx = contig ;
end ;
  
fname=genome_info.flat_fnames{contig_idx} ;
if genome_info.alphabet(1)=='a',
  NN='n' ;
else
  NN='N' ;
end ;

d=dir(fname) ; left_n = 0 ; right_n=0 ;
if isinf(stop), stop=d.bytes; end ;

if start<1 | stop>d.bytes,
  warning('load_genomic:contig_boundary', 'boundary of contig reached (start: %i, stop: %i, d.bytes: %i), padding with "n"', start, stop,d.bytes) ;

  if start<1,
    left_n = -(start-1) ;
    start = 1 ;
  end ;
  if stop>d.bytes,
    right_n = stop-d.bytes ;
  end ;
  if start>d.bytes && stop>d.bytes,
    str = char(NN*ones(1,stop-start)) ;
    return ;
  end ;
end ;

fd=fopen(fname,'r') ;
fseek(fd, start-1, -1) ;
str=char(fread(fd,stop-start+1, 'char=>char'))' ;
fclose(fd) ;
if left_n>0,
  str=char([NN*ones(1,left_n) str]) ;
end ;
if right_n>0,
  str=char([str NN*ones(1,right_n)]) ;
end ;

assert(length(str)==stop-start+left_n+1)
%let = bigunique(str) ;
let = unique(str) ;
if genome_info.alphabet(1)=='a' && ~isempty(setdiff(let,'acgtn')),
  % warning(sprintf('invalid letter in genome! ("%s")', setdiff(let,'acgtn')));
  let = setdiff(let, 'acgtn') ;
  for c=let,
    str(str==c) = 'n';
  end ;
elseif genome_info.alphabet(1)=='A' && ~isempty(setdiff(let,'ACGTN')),
  warning(sprintf('invalid letter in genome! ("%s")', setdiff(let,'ACGTN')));
  let = setdiff(let, 'ACGTN') ;
  for c=let,
    str(str==c) = 'N';
  end ;
end

if strand=='-',
  str=reverse_complement(str) ;
end ;

return ;
