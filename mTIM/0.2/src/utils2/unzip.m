function f=unzip(filename, dsp) 
% f=unzip(filename) 

if nargin<2,
  dsp=0 ;
end 
f=0 ;
if ~fexist(filename) & fexist([filename '.gz']),
  if dsp,
    disp(['ungzipping ' filename '.gz']) ;
  end ;
  unix(['gzip -d ' filename '.gz']) ;
  f=1 ;
end ;
