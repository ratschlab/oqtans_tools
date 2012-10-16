function str = splitstr(str,splt,ins);

% SPLITSTR Split String
%
% SPLITSTR( String , Split , Insert )
%
% Split = RowLength + i*MaxLength
%
% RowLength  Maximum Length of Row,
% MaxLength  Maximum Length of String
%
% RowLength < 0  ==> call of RMBLANK
% 
% RowLength + Tolerance/RowLength with  Tolerance <= RowLength
%
% Insert     Strings, after each a NewLine will insertet 
%

Nin = nargin;

nl = char(10);                         % NewLine  (InsertCharacter)
bl = { char(9)  char(32)  char(160) }; % Blanks
sc = { ' - '  '- '  '. '  '; '  ', '   ': '  ' ' }; % SplitCharacter

qt = '"';
ap = ' ... ';  % AppendString

%-----------------------------------------
% Check Inputs

if isempty(str)
   return
end

if Nin < 2
   splt = [];
end

if Nin < 3
   ins = [];
end

if ~chkstr(str)
    error('Input must be a String.');
end

if isempty(splt)
   splt = 60;
end

nmax = abs(imag(splt(1)));
splt =     real(splt(1));

if ~isempty(ins)
    [ok,ins] = chkcstr(ins);
    if ~ok
        ins = [];
    end
end

%-----------------------------------------
% Get NMAX and OVerSize from SPLiT

nmax = abs(imag(splt(1)));
splt =     real(splt(1));

rmbl = ( splt > 0 );

splt = abs(splt);

ovs  = splt;
splt = floor(splt);
ovs  = ovs - splt;
ovs  = 1e-10 * round(ovs/1e-10);
ovs  = splt + ceil( ovs * splt );

%-----------------------------------------
% Check with Insert

if ~isempty(ins)
    for cc = ins(:)'
        str = strrep(str,cc{1},[cc{1} nl]);
    end
end

%-----------------------------------------
% Check for Quota

n  = size(str,2);
iq = all( str([1 n]) == qt );
if iq
   str = str( 2 : n-1 );
    n  = n - 2;
else
   qt = '';
end

%-----------------------------------------
% RMBLANK if positive SPLIT

if rmbl
   str = rmblank(str,2);
end

if isempty(str)
   str = fill_str(str,[],qt);
   return
end

%-----------------------------------------
% Surround with NL

str = cat( 2 , nl , str , nl );

%-----------------------------------------
% Check for Blanks

[str,ib,lb] = chk_blk(str,bl,nl);

%-----------------------------------------
% Check with NMAX of String

n  = size(str,2);
na = size(ap,2);

nmax = nmax + n * ( nmax == 0 );
nmax = max( nmax , na+1 );

nmax = nmax + ( sum(str) == nl );

if all( n <= [ splt+n*(splt==0) nmax ]  );
   str = fill_str(str(2:n-1),bl);
   return
end

if n <= nmax

   ap = '';

else

   nmax = nmax - na;

   str = rmblank( str(1:nmax) , 2-i );

   if isempty(str)
      str = fill_str(str,[],qt,ap);
      return
   end

   n   = size(str,2);

   if ~isempty(ib)
       jj  = ( ib+lb-1 > n );
       if any(jj)
              jj  = find(jj);
           ib(jj) = [];
           lb(jj) = [];
       end
   end

end

%-----------------------------------------
% Check for SPLIT and 
% Length of Segments between NewLine

ret = ( splt == 0 );

if ~ret
    ret = cat( 2 , 0 , find(str==nl) , n+1 );
    ret = diff( ret , 1 , 2 ) - 1;   % Length between NL
    ret = all( ret <= splt );
end

if ret
   str = fill_str(str(2:n-1),bl,qt,ap);
   return
end

%-----------------------------------------
% Check for Characters to Split

[im,fm,jm] = get_mark( str , [{nl} sc] , 0 );

if all( fm == 1 )                % NL only !!!
   str = fill_str(str,bl,qt,ap);
   return
end

fm = fm - 1;    % NL == ZERO !!!

jm = im + jm(fm+1) - 1;  % EndIndex 

%-----------------------------------------
% 1: Remove Marker with same EndIndex
% 0: Remove Marker with following NL

for ii = [ 1  0 ]

    jj = ( 1 : size(im,2)-1 );

    if ii
       jj = ( ( jm(jj+1) == jm(jj) ) | ( im(jj+1) == im(jj) ) );
    else
       jj = ( ( fm(jj+1) == 0 ) & ~( fm(jj) == 0 ) & ...
              ( im(jj+1)-jm(jj) == 1 ) );
    end

    if any(jj)
           jj  = find(jj) + ii;
        im(jj) = [];
        jm(jj) = [];
        fm(jj) = [];
    end

end

%-----------------------------------------
% Set Length of NL-following Blanks

if isempty(ib);
   
   blk = '';
   ibl = [];
    jb = zeros(size(im));
    lb = jb;

else

     mb = size(ib,2);

    blk = str(grp2ind(ib+1,lb));        %  NL-following Blanks to insert
    ibl = cumsum(cat(2,1,lb(1:mb-1)));  % Index  in BLK

    nb     = zeros(1,n);
    nb(ib) = lb + i*( 1 : mb );

    ib = im(find(fm==0));
    lb = nb(ib);

    jb = imag(lb);  % Index  in IBL and NBL
    lb = real(lb);  % Length in BLK

    %----------------------------------------------------
    % Shift if previous Larger and NOT surrounded by ZEROS

    mb      = size(ib,2);

    ii     = zeros(1,mb);
       jj  = ( 2 : mb-1 );
    ii(jj) = ( ( lb(jj-1) == 0 ) & ( lb(jj+1) == 0 ) ); % Surrounded  by ZEROS

       jj  = ( 2 : mb );
    ii(jj) = ( ( lb(jj) < lb(jj-1) ) & ~ii(jj-1) );     % Shift if previous larger

    if any(ii)
          ii  = find(ii);
       jb(ii) = jb(ii-1);
       lb(ii) = lb(ii-1);
    end

    %----------------------------------------------------

    nb     = zeros(1,n);

       jj  = ( mb : -1 : 1 );

       ib  = n - ib(jj) + 1;  % Start from End
    nb(ib) = lb(jj) + i*jb(jj);

    jj = ( 1 : mb-1 );

    nb(ib(jj+1)) = nb(ib(jj+1)) - nb(ib(jj));
    nb = cumsum(nb,2);
    nb = nb( n : -1 : 1 );

    lb = nb(im);

    jb = imag(lb);  % Index  in IBL and NBL
    lb = real(lb);  % Length in BLK

end
 
%-----------------------------------------

fb = size(sc,2) - 1;             % FM of Blank

fm = ( fm > 0 ) + ( fm == fb );  % 0: NL; 2: BL; 1: other

nm  = size(im,2);

lm = diff(cat(2,jm(1),jm),1,2);    % Length between Characters to Split, incl.

lm = cumsum(lm,2);

ind = ( 2 : nm );

%%% [ splt  ovs ]

while 1

   % Remove Length-Offset

   ln     = zeros(1,nm);
      ii  = find( fm <= 0 );
   ln(ii) = 1;
   ln     = cumsum(ln,2);
   ln     = lm(ii(ln));

   ls = ( lm - ln ) .* ~isnan(fm);

   ns = floor(ls/splt);

   jj = ( ( ns(ind-1) == 0 ) & ( ns(ind) > 0 ) );

   if ~any(jj)
       break
   end

   jj = find(jj) + 1;

   % NewLine in Range of OverSize or Split after Blank in Range of OverSize
   kk = ( ( ( fm(jj+1) == 0 ) & ( ls(jj+1) <= ovs ) ) | ...
          ( ( fm(jj+1) == 1 ) & ( ls(jj+1) <= ovs ) & ( fm(jj) == 2 ) ) );

   %%% lk = 0*ok; lk(jj) = 1-2*kk; [i1 ; l1 ; ok ; n1 ; ls ; ns ; lk] 

   fm(jj) = -1;

   if any(kk)
        nm   = nm - sum(kk);
        ind  = ind( 1 : (nm-1) );
         kk  = jj(find(kk));
      jm(kk) = [];
      lm(kk) = [];
      fm(kk) = [];
      jb(kk) = [];
      lb(kk) = [];
   end

end

ok = ( fm == -1 );

if any(ok)

   ok = find(ok);

   jm = jm(ok);
   jb = jb(ok);
   lb = lb(ok);

   [str,jj] = insert(str,jm,nl);

   ok = ( ( lb > 0 ) & ( jb > 0 ) );

   if any(ok)
      ni = sum(ok);
      ok = find(ok);
      jj = jj(ok);
      jb = jb(ok);
      lb = lb(ok);
      for ii = ( ni : -1 : 1 )
          kk = ibl(jb(ii)) + ( 1 : lb(ii) ) - 1;
          str = cat( 2 , str(1:jj(ii)) , blk(kk) , str(jj(ii)+1:end) );
      end
   end

end

str = fill_str(str(2:end-1),bl,qt,ap);


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = fill_str(str,bl,qt,ap);

Nin = nargin;

if Nin < 2
   bl = [];
end
if Nin < 3
   qt = '';
end
if Nin < 4
   ap = '';
end

if ~isempty(bl)
    for ii = 1 : size(bl,2)
        str(find(str==ii)) = bl{ii};
    end
end

if ~isempty(ap)
    str = cat(2,str,ap);
end

if ~isempty(qt)
    str = cat(2,qt,str,qt);
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [str,ib,lb] = chk_blk(str,bl,nl);

n = size(str,2);

%-----------------------------------------
% Groups of Blanks

ib = zeros(1,n);
for cc = bl
    ib = ( ib | ( str == cc{1} ) );
end

[ib,lb] = ind2grp( find(ib) );

%-----------------------------------------
% Check for following NL

jj = ( str(ib+lb) == nl );  

if any(jj)

   jj = find(jj);

   ok = ~( str(ib(jj)-1) == nl );  % Previous is NOT NL

   ok = ok(:);

  lb(jj) = lb(jj) - ok;

  s( grp2ind(ib(jj)+ok,lb(jj)) ) = [];

   nb     = 0 * lb;
   nb(jj) = lb(jj);
   nb     = cumsum( cat(1,0,nb) , 1 );
   ib     = ib - nb(1:(end-1));

   if any(~ok)
            ok   = find(~ok);
      ib(jj(ok)) = [];
      lb(jj(ok)) = [];
   end

end

%-----------------------------------------
% Check for previous NL, Mask Blanks

jj = ( str(ib-1) == nl );  

if any(jj)

   jj = find(jj);

   ib = ib(jj);   % Index of first Blank after NL 
   lb = lb(jj);   % Number of following Blanks

   jj = grp2ind(ib,lb);

   for ii = 1 : size(bl,2)
       kk = find( str(jj) == bl{ii} );
       str(jj(kk)) = ii;
   end

   ib = ib - 1;   % Index of NL

   ib = ib(:)';
   lb = lb(:)';

else

   ib = zeros(1,0);
   lb = zeros(1,0);

end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [im,fm,lm] = get_mark(str,mark,ret);

% GET_MARK  Returns MarkerIndize, Flag and Length
%

nm = size(mark,2);

im = zeros( 1 , 0  );  % Index for Marker
fm = zeros( 1 , 0  );  % Flag  for Marker
lm = zeros( 1 , nm );  % Length of Marker

%---------------------------------------------------------------------
% Find Marker: Index & Flag

for ii = 1 : nm

    lm(ii) = size( mark{ii} , 2 );

    if ~( lm(ii) == 0 )

        jj = findstr( str , mark{ii} );

        if ~isempty(jj)

            im = cat( 2 , im ,      jj );
            fm = cat( 2 , fm , ii+0*jj );
    
        elseif any( ii == ret )

            return
 
        end

    end

end

%---------------------------------------------------------------------
% Sort

[im,si] = sort(im);

fm = fm(si);


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = rmblank(str,dim,cc)

% RMBLANK  Remove Blanks, NewLines at Begin and End of CharacterArrays
%
% String = RMBLANK( CharArray )
%
% CharArray  2-dimensional CharacterArray
%
% further Options:
%
% String = RMBLANK( CharArray , DIM , CHAR )
%
%  
%  DIM  specifies Dimension to work, 
%       default: 2
%
%    A positive complex Part of DIM, to removes Blanks only from Start,
%    A negative complex Part of DIM, to removes Blanks only from End.
%       
%  CHAR specifies BlankCharacters to remove
%       default:  [ 160  32  13  10  9  0 ];  % [ NBSP Space CR LF TAB ZERO ]
%

  
msg = '';
 nl = char(10);

Nin = nargin;

if Nin < 1
  error('Not enough Input Arguments.')
else
  if ischar(str)
    str = double(str);
  end
  ok = isnumeric(str);
  if ok
    ok = all( ( mod(str(:),1) == 0 )  & ...
              ( str(:) >= 0 ) & isfinite(str(:))  );
  end
  if ~ok
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CharArray must be a String or ASCII-Codes.'];
  end
  if size(str,1)*size(str,2) ~= prod(size(str))
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CharArray must be 2-dimensional.'];
  end     
end

if Nin < 2
  dim = 2;
else
  if ~isnumeric(dim)
    msg = [ msg nl(1:(end*(~isempty(msg)))) ...
            'Input DIM must be numeric.' ];
  elseif ~isempty(dim)
    dim = dim(:);
    if ~all( ( real(dim) == 1 ) |  ( real(dim) == 2 ) )
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
             'Values for Input DIM must define 1. or 2. Dimension.' ];
    end
  end 
end

if Nin < 3
  cc = [ 160  32  13  10  9  0 ];  % [ NBSP  Space CR LF TAB ZERO ]
else
  if ischar(cc)
    cc = double(cc);
  end
  ok = isnumeric(cc);
  if ok & ~isempty(cc)
    cc = cc(:)';
    ok = all( ( mod(cc,1) == 0 )  & ...
              ( cc >= 0 ) & isfinite(cc)  );
  end
  if ~ok
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CHAR must be a String or ASCII-Codes.'];
  end
end

if ~isempty(msg)
  error(msg)
end


if isempty(str)
 str = '';
 return
end

if isempty(dim) | isempty(cc)
  str = double(str);
  return
end


  blank  = 0*str;

  for ii = cc
    blank = ( blank | ( str == ii ) );
  end

  si = size(str);

  for ii = 1 : size(dim,1)

    d = dim(ii);

    s = sign(imag(d));  % Remove from wich Side:  1  0  -1 
 
    d = real(d);

    jj = find( sum(blank,3-d) == si(3-d) );  % Columns with full Blanks

    if ~isempty(jj) 

         p  = [ 3-d  d ];
        str = permute(str,p);

         jj = jj(:)';
         nb = size(jj,2);

        %--------------------------------------------
        % Blank at Begin

        ind = ( 1 : nb );
        jj1 = find( ( ( jj == ind ) & ( s >= 0 ) ) );

        %--------------------------------------------
        % Blank at End

        ind = ind + si(d) - nb;
        jj2 = find( ( ( jj == ind ) & ( s <= 0 ) ) );

        %--------------------------------------------

        str(:,jj([jj1 jj2])) = [];

        str = permute(str,p);

    end
    
  end

  str = char(str);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [i0,l] = ind2grp(ii);

% IND2GRP  Built StartIndex and Length from IndexVector
%
% [ StartIndex , GroupLength ] = IND2GRP( Index )
%

i0 = zeros(0,1);
l  = zeros(0,1);

if isempty(ii);
   return
end

ii = ii(:);
n  = size(ii,1);

i0 = cat( 1 , 1 , find( diff(ii,1,1) > 1 )+1 , n+1 );

l  = diff(i0,1,1);

i0 = ii(i0(1:end-1));


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ii = grp2ind(i0,l,s);

% GRP2IND  Built IndexVector from StartIndex and Length
%
% Index = GRP2IND( StartIndex , GroupLength )
%
% Index = GRP2IND( StartIndex , GroupLength , LowSampleStep  )
%

ii = [];

if isempty(i0);
   return
end

if nargin < 3
   s = 1;
end

si = size(i0);

if ( sum( si > 1 ) > 1 )
   error('StartIndex must be a Vector.');
end

i0 = i0(:);
l  =  l(:);

if ~isequal(size(i0,1),size(l,1))
   error('Size of StartIndex and GroupLenght must be the same.');
end

jj = ( l == 0 );
if all(jj)
   return
end

if any(jj)
      jj  = find(jj);
   i0(jj) = [];
    l(jj) = [];
end

n = size(l,1);

l = ceil( l / s );

ii = s * ones(sum(l),1);
jj = cumsum( cat(1,1,l) , 1 );

ii(jj(1:n)) = i0;

if n > 1
   ii(jj(2:n)) = ii(jj(2:n))-(i0(1:n-1,1)+s*(l(1:n-1)-1));
end

ii = cumsum(ii,1);

if n == 1
   return
end

jj = find( si > 1 );

if jj == 1
   return
end
  
perm = cat( 2 , (1:jj-1)+1 , 1 , (jj+1:size(si,2)) );

ii = permute(ii,perm);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkstr(str,opt)


% CHKSTR  Checks Input for String
%
%  ok = chkstr(str,Option)
%
%  Option ~= 0 ==> only true for nonempty Strings
%
%   default:   Option == 0  (true for empty Strings)
%

 
if nargin < 2
   opt = 0;
end

ok = ( strcmp( class(str) , 'char' )      & ...
       ( prod(size(str)) == size(str,2) ) & ...
       ( isequal(opt,0) | ~isempty(str) )         );


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,str] = chkcstr(str,opt)


% CHKCSTR  Checks Input for CellString, contains Strings !
%
%  [ok,str] = chkcstr(str,Option)
%
%  Option ~= 0 ==> CharacterArrays not allowed,
%
%   default: Option == 0   ==>  CharacterArrays --> CellString
%
 
if nargin < 2
   opt = 0;
end

if strcmp(class(str),'char') & isequal(opt,0)
   n = size(str,1);
   if n == 1
      str = strrep(str,char(32),char(1));
   end
   str = cellstr(str);
   if n == 1
      str = strrep(str,char(1),char(32));
   end
end

ok = iscellstr(str);
if ~ok
   return
end

try
  s = cat(2,str{:});
catch
  ok = 0;
  return
end
 
ok = ( strcmp(class(s),'char')  &  ( prod(size(s)) == size(s,2) ) );

