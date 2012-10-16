function k= kurtosisP(signalMatrix)
% function k= kurtosisP(signal)

kk=size (signalMatrix) ;

numSignals = min(kk);
numData=max(kk);

if kk(1)>kk(2),
  signalMatrix=signalMatrix';
end

signalMatrix = signalMatrix - mean(signalMatrix')' * ones(1,numData);
signalMatrix = signalMatrix ./ (std(signalMatrix')' * ones(1,numData));
  
k = mean((signalMatrix.^4)') - 3 * (mean((signalMatrix.^2)').^2);
