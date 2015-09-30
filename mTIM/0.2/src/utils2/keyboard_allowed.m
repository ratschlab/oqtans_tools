function ret=keyboard_allowed()
% ret=keyboard_allowed()
%
% returns 1, if a keyboard command would be allowed

global g_ignore_keyboard
global THIS_IS_A_RPROC_PROCESS  

if isequal(g_ignore_keyboard, 1),
  ret=0 ;
  return ;
end ;

if isequal(THIS_IS_A_RPROC_PROCESS, 1),
  ret=0 ;
  return ;
end ;

ret=1 ;

