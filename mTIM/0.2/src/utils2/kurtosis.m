function kurt=kurtosis(data)

[m,n]=size(data);

if min(size(data))==1,
data=data(:);
end

%data=data-ones(m,1)*mean(data);

kurt=mean(data.^4)./(diag(cov(data)).^2)'-3;
