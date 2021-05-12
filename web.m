function [ranges] = web(singleEps,peaks,indices,step,ax,dat0,param,sigma,limit,func,maxIter)

len = length(indices);
ranges = [singleEps singleEps];
if maxIter == 0
    return
end

for i = 1:len
    mask = zeros(6,1);
    mask(indices(i)) = 1;
    
    manyEpsTest = singleEps + [-1 1] * step .* mask;
    chi2Test = func(manyEpsTest,dat0,param,peaks)/(sigma^2);
    
    if any(chi2Test < limit)
    
        manyEps = singleEps + ax .* mask;

        chi2 = func(manyEps,dat0,param,peaks)/(sigma^2);

        smallPoints = chi2 < limit;
        findFirst = find(smallPoints, 1, 'first');
        findLast = find(smallPoints, 1, 'last');

        ranges(indices(i),1) = min(ranges(indices(i),1),manyEps(indices(i),max(findFirst - 1,1)));
        ranges(indices(i),2) = max(ranges(indices(i),2),manyEps(indices(i),min(findLast + 1,length(ax))));
        
        newIndices = 1:6;
        newIndices(indices(i)) = [];
        
        firstEps = singleEps;
        firstEps(indices(i)) = manyEps(indices(i),findFirst);
        firstRanges = web(firstEps,peaks,newIndices,step,ax,dat0,param,sigma,limit,func,maxIter - 1);
        ranges(:,1) = min(ranges(:,1), firstRanges(:,1));
        ranges(:,2) = max(ranges(:,2), firstRanges(:,2));
        
        if findFirst ~= findLast
            lastEps = singleEps;
            lastEps(indices(i)) = manyEps(indices(i),findLast);
            lastRanges = web(lastEps,peaks,newIndices,step,ax,dat0,param,sigma,limit,func,maxIter - 1);
            ranges(:,1) = min(ranges(:,1), lastRanges(:,1));
            ranges(:,2) = max(ranges(:,2), lastRanges(:,2));
        end
    end
end

end