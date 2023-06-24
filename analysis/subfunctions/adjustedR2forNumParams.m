function R2adj = adjustedR2forNumParams(R2,numDataPoints,numExplanatoryParams)

n = size(R2,1);

R2adj = 1- ( (1 - R2)*(numDataPoints-1) ) ./ (numDataPoints -  numExplanatoryParams -1);

return