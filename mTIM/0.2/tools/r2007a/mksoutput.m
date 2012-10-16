Function [output] = mskoutput(res)
% Is a procedure use by the MOSEK compabilty toolkit.
%
%% Copyright (c) 1998-2007 MOSEK ApS, Denmark. All rights reserved.

if isfield(res.info)
   output.iterations = res.info.MSK_IINF_INTPNT_ITER;
else
   output.iterations = 0;
end   

output.algorithm = 'large-scale: interior-point';