function list=who_file(fname) ;
% list=who_file(fname) ;
  
engine = determine_engine() ;
if isequal(engine, 'matlab')
  list=who('-file', fname) ;
else
  L=load(fname) ;
  list=fieldnames(L) ;
end ;
