function compare_structures(a,b, verb) 
% compare_structures(a,b) 

af=fields(a) ;
bf=fields(b) ;
if nargin<3, verb=0 ; end ;

da=setdiff(af,bf) ;
if ~isempty(da),
  fprintf('fields only in a:\n') ;
  da 
end ;
db=setdiff(bf,af) ;
if ~isempty(db),
  fprintf('fields only in b:\n') ;
  db 
end ;

f=intersect(af,bf);
for i=1:length(f),
  if ~isequalwithequalnans(a.(f{i}), b.(f{i})),
    fprintf('field %s differs:\n', f{i})
    if verb>0,
      a.(f{i})
      b.(f{i})
    end ;
  end 
end 
