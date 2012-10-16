function smRowVectors = softmax(RowVectors)
% function smRowVectors = softmax(RowVectors)

%RowVectors(:,1:10) 
eRV=exp(RowVectors) ;
%maxmax=max(max(eRV))
Di=ones(size(eRV)) ;
i=0 ;
for eV=eRV
	i=i+1 ;
	Di(:,i)=Di(:,i)*(1+sum(eV)) ;
end ;
smRowVectors=eRV./Di ;

%smRowVectors(:,1:10) 
   