function f=zip(filename, dsp, bg) 
%f=zip(filename) 

if nargin<2,
  bg=0 ;
end ;

f=0 ;
if fexist(filename) & ~fexist([filename '.gz']),
  if bg,
    if dsp,
      disp(['gzipping ' filename '.gz in background']) ;
    end ;
    unix(['gzip ' filename ' &']) ;
  else
    if dsp,
      disp(['gzipping ' filename '.gz']) ;
    end ;
    unix(['gzip ' filename ]) ;
  end ;
  f=1 ;
end ;
