function [r,res] = mosekopt(cmd,prob,param,callback)
%
% Syntax     : [r,res] = mosekopt(cmd,prob,param,callback)
%
%
% Purpose    : Interface to the MOSEK optimization tools.
%          
%
% Description: Required arguments.
%                cmd       Commands e.g. 'minimize'.
%                prob      The optimization problem.   
%
%              Optional arguments.
%                param     MOSEK parameters.
%                callback  A callback function.
%
%
%              Please see "The MOSEK optimization toolbox manual"
%              for a detailed description of mosekopt. 
%
%
% Notes      : 1) Log information can be turned off by including echo(0)
%                 in cmd e.g mosekopt('minimize echo(0)',prob).
%
%              2) This is only a help file.   
%
%% Copyright (c) 1998-2007 MOSEK ApS, Denmark. All rights reserved.


