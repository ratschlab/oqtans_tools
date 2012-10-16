function mskstatus(callproc,verb,dual,rcode,res)
% Internal function used by linprog, quadprog, etc.
%
%% Copyright (c) 1998-2007 MOSEK ApS, Denmark. All rights reserved.

if ( isfield(res,'symbcon') )
    sc = res.symbcon;
else    
    [r,res] = mosekopt('symbcon');
    sc      = res.symbcon;
end

switch ( rcode )
case { sc.MSK_RES_ERR_INV_PROBLEM }
   disp([callproc ': Invalid problem which MOSEK cannot handle.']);
   disp([callproc ': Probably the problem is nonconvex.']);
end

if verb > 0
   if ( isfield(res,'sol') )
      if ( ( ~dual & res.sol.itr.prosta==sc.MSK_PRO_STA_PRIM_INFEAS ) | ...
           ( dual & res.sol.itr.prosta==sc.MSK_PRO_STA_PRIM_AND_DUAL_INFEAS ) )
         disp('The primal problem is infeasible.');
      end
      if ( ( ~dual & res.sol.itr.prosta==sc.MSK_PRO_STA_DUAL_INFEAS ) | ...
           ( dual & res.sol.itr.prosta==sc.MSK_PRO_STA_PRIM_INFEAS ) )
         disp('The dual problem is infeasible.');
      end
   end   
   
   if ( mskeflag(rcode,res)==1 )
      disp('Optimization terminated successfully.');   
   end
   
   if ( mskeflag(rcode,res)==0 )
     disp('Maximum number of iterations exceeded.')
   end
end
