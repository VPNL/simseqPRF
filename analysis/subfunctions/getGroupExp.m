
function [groupExp, groupMaxExpVal,groupMaxExpIdx] = getGroupExp(x,y)
        binsExp = 0.1:0.1:1;
     
        [groupMaxExpVal,groupMaxExpIdx] = max( reshape(y(:,x,:), size(y,1),[]));
        groupExp = binsExp(groupMaxExpIdx);
end