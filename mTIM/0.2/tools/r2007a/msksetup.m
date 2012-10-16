function [cmd,verbosity,param] = msksetup(domin,options)
% Purpose: Function used by the MOSEK compability toolbox.
%
%% Copyright (c) 1998-2007 MOSEK ApS, Denmark. All rights reserved. 

verbosity = 0;

param     = [];  

echolev   = 0;
prlev     = 0;                        

if ~isempty(options) & isa(options,'struct') & isfield(options,'Display')
  dispstr = deblank(getfield(options,'Display')); 
  if strcmp(dispstr,'off')
    echolev                           = 0;
  elseif strcmp(dispstr,'final') 
    echolev                           = 3;
    verbosity                         = 1;
  elseif strcmp(dispstr,'iter')  
    param.MSK_IPAR_LOG_INTPNT         = 1;
    param.MSK_IPAR_LOG_SIM            = 1;
    echolev                           = 3;
    verbosity                         = 1;
  else  
    echolev                           = 0;
  end  
end   
 
if isa(options,'struct') & isfield(options,'Diagnostics')
  diagstr = lower(getfield(options,'Diagnostics'));
  if strcmp(diagstr,'on') 
     fprintf('Diagnostics is on\n'); 
     echolev   = max(echolev,10);      
     prlev     = max(prlev,1); 
     verbosity = 1;
  else
     param.MSK_IPAR_MAX_NUM_WARNINGS = 0;
  end   
else
   param.MSK_IPAR_MAX_NUM_WARNINGS = 0;
end  

if echolev==0
   param.MSK_IPAR_LOG = 0;
end  
param.MSK_IPAR_LOG_INTPNT    = prlev;
param.MSK_IPAR_LOG_SIM       = prlev;
param.MSK_IPAR_LOG_BI        = prlev;
param.MSK_IPAR_LOG_PRESOLVE  = prlev;

%switch ( optimget(options,'diagnostics') )
%  case 'on' 
%     echolev = max(echolev,1);      
%  otherwise
%     param.MSK_IPAR_MAX_NUM_WARNINGS = 0;
%end % switch

% fprintf('@ done\n');

%
% Copy Mosek options.
%

optnames = fieldnames(options);
[m,n]    = size(optnames);

for i=1:m
  if ( strncmp('MSK_',optnames(i),4) )
    param = setfield(param,optnames{i},getfield(options,optnames{i}));    
  end   
end % for  

if ( domin )
   cmd = sprintf('minimize info echo(%d) statuskeys(1) symbcon',echolev);        
else
   cmd = sprintf('maximize info echo(%d) statuskeys(1) symbcon',echolev);
end
