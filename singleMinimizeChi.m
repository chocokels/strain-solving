function [sol, chi, peaksMeas, flag] = singleMinimizeChi(peaks,peaksErr,startPoints,dat0,param)

% For reference:
% strain = [exx,eyy,ezz,exy,eyz,ezx];

numPts = size(peaks,3);

exx = @(e)[ ...
    1/3 * e(1,:) - 2*sqrt(2)/3 * e(6,:) + 2/3 * e(3,:);
    1/3 * e(1,:) + 2*sqrt(2)/3 * e(6,:) + 2/3 * e(3,:);
    1/3 * e(2,:) + 2*sqrt(2)/3 * e(5,:) + 2/3 * e(3,:);
    1/3 * e(2,:) - 2*sqrt(2)/3 * e(5,:) + 2/3 * e(3,:)];
eyy = @(e)[ ...
    e(2,:);
    e(2,:);
    e(1,:);
    e(1,:)];
ezz = @(e)[ ...
    2/3 * e(1,:) + 2*sqrt(2)/3 * e(6,:) + 1/3 * e(3,:);
    2/3 * e(1,:) - 2*sqrt(2)/3 * e(6,:) + 1/3 * e(3,:);
    2/3 * e(2,:) - 2*sqrt(2)/3 * e(5,:) + 1/3 * e(3,:);
    2/3 * e(2,:) + 2*sqrt(2)/3 * e(5,:) + 1/3 * e(3,:)];
exy = @(e)[ ...
    sqrt(3)/3 * e(4,:) - sqrt(6)/3 * e(5,:);
    sqrt(3)/3 * e(4,:) + sqrt(6)/3 * e(5,:);
    - sqrt(3)/3 * e(4,:) - sqrt(6)/3 * e(6,:);
    - sqrt(3)/3 * e(4,:) + sqrt(6)/3 * e(6,:)];
eyz = @(e)[ ...
    sqrt(6)/3 * e(4,:) + sqrt(3)/3 * e(5,:);
    sqrt(6)/3 * e(4,:) - sqrt(3)/3 * e(5,:);
    - sqrt(6)/3 * e(4,:) + sqrt(3)/3 * e(6,:);
    - sqrt(6)/3 * e(4,:) - sqrt(3)/3 * e(6,:)];
ezx = @(e)[ ...
    sqrt(2)/3 * e(1,:) - 1/3 * e(6,:) - sqrt(2)/3 * e(3,:);
    sqrt(2)/3 * e(1,:) + 1/3 * e(6,:) - sqrt(2)/3 * e(3,:);
    sqrt(2)/3 * e(2,:) + 1/3 * e(5,:) - sqrt(2)/3 * e(3,:);
    sqrt(2)/3 * e(2,:) - 1/3 * e(5,:) - sqrt(2)/3 * e(3,:)];

eZPL = @(e) dat0.ZPL0 ...
    + param.tpara * ezz(e) + param.tperp * (exx(e) + eyy(e));
eDeltaEs = @(e) sqrt(dat0.deltaEs0^2 ...
    + 4 * (param.des * (exx(e) - eyy(e)) + param.fes * ezx(e)).^2 ...
    + 4 * (-2 * param.des * exy(e) + param.fes * eyz(e)).^2);
eDeltaGs = @(e) sqrt(dat0.deltaGs0^2 ...
    + 4 * (param.dgs * (exx(e) - eyy(e)) + param.fgs * ezx(e)).^2 ...
    + 4 * (-2 * param.dgs * exy(e) + param.fgs * eyz(e)).^2);

ePeaks = @(e) permute(eZPL(e),[1 3 2]) + permute(eDeltaEs(e),[1 3 2]) .* [-1 -1 1 1] / 2 + permute(eDeltaGs(e),[1 3 2]) .* [-1 1 -1 1] / 2;

sumSq = @(e,i) sum((peaks(:,:,i) - ePeaks(e)).^2 ./ (peaksErr(:,:,i)).^2,[1 2]);
sol = zeros(6,numPts);
chi = zeros(1,numPts);
peaksMeas = zeros(4,4,numPts);

flag = ones(numPts,1);

for i = 1:numPts
    sumSqSub = @(e) sumSq(e,i);
    options = optimset('TolX',1e-4,'TolFun',1e-8); %'PlotFcns',@optimplotfval, %'Display','iter'
    [sol(:,i),~,flagSingle] = fminsearch(sumSqSub,startPoints(:,i),options);
    flag(i) = flagSingle;
    chi(i) = sumSqSub(sol(:,i));
    peaksMeas(:,:,i) = ePeaks(sol(:,i));
end

end