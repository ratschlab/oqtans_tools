function save_append(fname__, append_flag, varargin)
% save_append(fname__, append_flag, varargin)
%
% octave compatible version of save -append

[engine] = determine_engine() ;

if isequal(engine, 'matlab')
  if length(varargin)>1,
    if append_flag,
      assert(mod(length(varargin),2)==0) ;
      for i=1:2:length(varargin),
        eval(sprintf('L__.%s = varargin{i+1};',  varargin{i})) ;
      end ;
      if exist(fname__, 'file')||exist([fname__ '.mat'], 'file')
      	save(fname__, '-V7',  '-append', '-struct', 'L__') ;
      else
      	save(fname__, '-V7', '-struct', 'L__') ;
	  end
      return ;
    end ;
  end ;
end ;

try,
  if append_flag
    L__=load(fname__) ;
  else
    L__=struct ;
  end ;
catch 
  L__=struct ;
end 

if length(varargin)>1,
  
  % assume varargin to have the following structure: varname1, var1, ...
  
  assert(mod(length(varargin),2)==0) ;
  for i=1:2:length(varargin),
    L__.(varargin{i}) = varargin{i+1} ;
  end ;
  
  f=fieldnames(L__) ;
  for i=1:length(f),
    eval(sprintf('%s=L__.%s;', f{i}, f{i})) ;
  end ;
  
  eval(sprintf('save -V7 %s %s', fname__, sprintf('%s ', f{:}))) ;
else
  
  % assume the second argument is a struct
  assert(length(varargin)==1) ;
  
  f = fieldnames(varargin{1}) ;
  for i=1:length(f),
    L__.(f{i}) = varargin{1}.(f{i}) ;
  end ;
  
  f=fieldnames(L__) ;
  for i=1:length(f),
    eval(sprintf('%s=L__.%s;', f{i}, f{i})) ;
  end ;

  eval(sprintf('save -V7 %s %s', fname__, sprintf('%s ', f{:}))) ;
  
end ;



