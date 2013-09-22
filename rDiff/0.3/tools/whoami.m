function user_name = whoami
% WHOAMI   Determines user name.
%
%   user_name = whoami
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2010 Gunnar Raetsch
%   Copyright (C) 2009-2010 Max Planck Society
%


global whoami_user_name  
  
if isempty(whoami_user_name)
  
  [ret, whoami_user_name]=unix('whoami') ; 
  if ret~=0,
    if keyboard_allowed(),
      keyboard ;
    end ;
    whoami_user_name='unknown' ;
  end ;
end ;

user_name= whoami_user_name ;
idx=find(user_name<30, 1, 'first') ;
user_name = user_name(1:idx-1) ;
