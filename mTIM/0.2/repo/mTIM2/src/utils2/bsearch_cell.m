% bsearch_cell(x,var)
% Written by Aroh Barjatya, adapted by Gunnar Raetsch to cell
% string arrays
% Binary search for values specified in vector 'var' within data
% vector 'x'
% The data has to be pre-sorted in ascending or decending order
% There is no way to predict how the function will behave if there 
% are multiple numbers with same value.
% returns the index values of the searched numbers

function index = bsearch_cell(x, var)
% index = bsearch_cell(x, var)

xLen = length(x) ;

for i = 1:length(var)
    low = 1;
    high = xLen;
    if strcomp(var{i}, x{low})<1
        index(i) = low;
        continue;
    elseif strcomp(var{i}, x{high})>-1
        index(i) = high;
        continue;
    end
    flag = 0;
    while (low <= high)
        mid = round((low + high)/2);
        if strcomp(var{i}, x{mid})==-1,
            high = mid;
        elseif strcomp(var{i}, x{mid})==1,
            low = mid;
        else
            index(i) = mid;
            flag = 1;
            break;
        end
        if (low == high - 1)
            break
        end
    end
    if (flag == 1)
        continue;
    end
    if (low == high)
        index(i) = low;
%    elseif ((x(low) - var(i))^2 > (x(high) - var(i))^2)
%        index(i) = high;
    else
        index(i) = low;
    end
end


function ret=strcomp(a,b) ;

if isequal(a,b),
  ret = 0 ;
else
  [tmp,idx]=sort({a,b}) ;
  if idx(1)==1,
    ret=-1 ;
  else
    ret=1 ;
  end ;
end ;
