function un=uname
%un=uname
%
% global SYSTEM_UNAME

global SYSTEM_UNAME

if equal(class(SYSTEM_UNAME), 'char')
  un=SYSTEM_UNAME ;
  return ;
end ;

% get host name 
!uname > /tmp/uname.log
hnd=fopen('/tmp/uname.log') ;
un=fscanf(hnd, '%s') ;
fclose(hnd) ;
!rm /tmp/uname.log

SYSTEM_UNAME=un 