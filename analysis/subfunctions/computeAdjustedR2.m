function R2adj = computeAdjustedR2(data,prediction,numExplanatoryParams)

R2 = computeCoD(data,prediction);
n = size(data,1);

R2adj = ( (1 - R2)*(n-1) ) ./ (n -  numExplanatoryParams -1);