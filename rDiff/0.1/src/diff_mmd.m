function [pval,info] = diff_mmd(reads1,reads2,gene,kernel)
% simple application of mmd to test for differential distributions
% of reads1, reads2
% reads1: N1 x L
% reads2: N2 x L
% note: if N1!=N2 the larger one is subsampled to yield the same
% size. This is needed to apply mmd
% kernel: currently polynomial only

%1. ensure reads are sparse

%reads1temp = sparse(reads1(:,sum([reads1;reads2],1)>0));
%reads2 = sparse(reads2(:,sum([reads1;reads2],1)>0));
%reads1=reads1temp;
bootstraps = 1000;

if nargin<4
  kernel = 'poly';
end;

N1 = size(reads1,1);
N2 = size(reads2,1);
N = N1 + N2;

if not(or(N1==0,N2==0))
statistic = eucl_dist(mean(reads1,1),mean(reads2,1))^2;

allreads = [reads1;reads2];


%Use the transpose to make the selection of columms faster
all_reads_trans=allreads';
for i = 1:bootstraps

  r = randperm(N);
  
  sample1 = all_reads_trans(:,r(1:N1));
  sample2 = all_reads_trans(:,r(N1+1:N));
  
  bootstrap_results(i) = eucl_dist(mean(sample1,2), mean(sample2,2))^2;
   
end

pval =  double(sum(bootstrap_results > statistic)) / bootstraps;
info = [];

else
pval = 1;
info = [];
  
end



%
%
%%subsample from the larger one:
%N1 = size(reads1,1);
%N2 = size(reads2,1);
%if(N1<N2)
%  Ir2 = randperm(N2);
%  Ir2 = Ir2(1:N1);
%  reads2 = reads2(Ir2,:);
%elseif (N2<N1)
%  Ir1 = randperm(N1);
%  Ir1 = Ir1(1:N2);
%  reads1 = reads1(Ir1,:);
%end;
%
%x = reads1;
%y = reads2;
%global Kxx Kyy Kxy                %kernel matrices
%if strcmp(kernel,'poly')
%  %polynomial kernel
%  Kxx = (x*x' + 1).^2;
%  Kyy = (y*y' + 1).^2;
%  Kxy = (x*y' + 1).^2;
%elseif strcmp(kernel,'exp')
%  %squared exp. or Gauss kernel
%  G = x*x';
%  L = size(G,1);
%  nor = G(1:L+1:L^2);
%  Kxx = -2*G + repmat(nor',[1,L]) +  repmat(nor,[L,1]);
%  clear G L nor 
%  % Kyy
%  G = y*y';
%  L = size(G,1);
%  nor = G(1:L+1:L^2);
%  Kyy = -2*G + repmat(nor',[1,L]) +  repmat(nor,[L,1]);
%  clear G L nor 
%  % Kxy
%  G = x*y';
%  L = size(G,1);
%  norx = sum(x.*x,2);
%  nory = sum(y.*y,2);
%  Kxy = (-2*G + repmat(norx,[1,length(nory)]) +  repmat(nory',[length(norx),1]))';
%  clear G L norx nory 
%  %now get the median distance
%  mdist = median(Kxy(Kxy~=0));
%  sigma = sqrt(mdist/2);
%  if sigma ==0
%    sigma =1;
%  end
%  %apply RBF
%  Kxx = exp(-1/2/sigma^2 * Kxx);
%  Kyy = exp(-1/2/sigma^2 * Kyy);
%  Kxy = exp(-1/2/sigma^2 * Kxy);
%end;
%%generate labels and data and run MMD
%labels = [ones(size(x,1),1);-ones(size(y,1),1)];
%[H,info] = mmd([x;y],labels,0.05,'bootstrap');         
%clear global  Kxx Kyy Kxy 
%%calcualte p-value from permutations manually:
%value = info.mmd.bootstrap.val;
%bootval = info.mmd.bootstrap.info.bootval;
%pval = (sum(value<bootval)+1)/size(bootval,1);
%info.H = H;
%

function result = eucl_dist(A,B)

result = sqrt(sum( (A - B) .^ 2 ));
