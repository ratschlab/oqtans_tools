function [pval,info] = rDiff_mmd(CFG,reads1,reads2)
% simple application of mmd to test for differential distributions
% of reads1, reads2
% reads1: N1 x L
% reads2: N2 x L


bootstraps=CFG.bootstraps;


% ensure reads are sparse and remove zero collumns
%reads1temp = sparse(reads1(:,sum([reads1;reads2],1)>0));
%reads2 = sparse(reads2(:,sum([reads1;reads2],1)>0));
%reads1=reads1temp;

statistic = eucl_dist(mean(reads1,1),mean(reads2,1))^2;

allreads = [reads1;reads2];

N1 = size(reads1,1);
N2 = size(reads2,1);
N = N1 + N2;

 
%Use the transpose to make the selection of columms faster
all_reads_trans=allreads';

%bootstraping
for i = 1:bootstraps

  r = randperm(N);
  
  sample1 = all_reads_trans(:,r(1:N1));
  sample2 = all_reads_trans(:,r(N1+1:N));
  
  bootstrap_results(i) = eucl_dist(mean(sample1,2), mean(sample2,2))^2;
   
end

%Calculate the p-value
pval =  min(1,double(1+sum(bootstrap_results >= statistic)) / bootstraps);
info = [];



function result = eucl_dist(A,B)

result = sqrt(sum( (A - B) .^ 2 ));
