function emacs(p1, p2)
%
if nargin==0,
  !emacs &
elseif nargin==1,
  eval(['!emacs ' p1]) ;
else
  eval(['!emacs v' p1 ' ' p2]) ;
end ;