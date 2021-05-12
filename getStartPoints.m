function [startPoints,flag] = getStartPoints(peaks,posax,dat0,param)
%GETSTARTPOINTS First step of strain solving
%   Simplifies to two family case and solves mathematically

exx = @(e)[ ...
    1/3 * e(1) + 2/3 * e(3);
    1/3 * e(2) + 2/3 * e(3)];
eyy = @(e)[ ...
    e(2);
    e(1)];
ezz = @(e)[ ...
    2/3 * e(1,:) + 1/3 * e(3,:);
    2/3 * e(2,:) + 1/3 * e(3,:)];
exy = @(e)[ ...
    sqrt(3)/3 * e(4,:);
    - sqrt(3)/3 * e(4,:)];
eyz = @(e)[ ...
    sqrt(6)/3 * e(4,:);
    - sqrt(6)/3 * e(4,:)];
ezx = @(e)[ ...
    sqrt(2)/3 * e(1,:) - sqrt(2)/3 * e(3,:);
    sqrt(2)/3 * e(2,:) - sqrt(2)/3 * e(3,:)];

eZPL = @(e) dat0.ZPL0 ...
    + param.tpara * ezz(e) + param.tperp * (exx(e) + eyy(e));
eDeltaEs = @(e) sqrt(dat0.deltaEs0^2 ...
    + 4 * (param.des * (exx(e) - eyy(e)) + param.fes * ezx(e)).^2 ...
    + 4 * (-2 * param.des * exy(e) + param.fes * eyz(e)).^2);
eDeltaGs = @(e) sqrt(dat0.deltaGs0^2 ...
    + 4 * (param.dgs * (exx(e) - eyy(e)) + param.fgs * ezx(e)).^2 ...
    + 4 * (-2 * param.dgs * exy(e) + param.fgs * eyz(e)).^2);

dat.ZPL = squeeze(mean(peaks,2));
dat.deltaEs = squeeze(mean(peaks(:,[3 4],:),2) - mean(peaks(:,[1 2],:),2));
dat.deltaGs = squeeze(mean(peaks(:,[2 4],:),2) - mean(peaks(:,[1 3],:),2));

dat2.ZPL = [mean(dat.ZPL(1:2,:)); mean(dat.ZPL(3:4,:))];
dat2.deltaEs = [mean(dat.deltaEs(1:2,:)); mean(dat.deltaEs(3:4,:))];
dat2.deltaGs = [mean(dat.deltaEs(1:2,:)); mean(dat.deltaEs(3:4,:))];

solEpsEs = zeros(4,length(posax));
solEpsGs = zeros(4,length(posax));
eps = sym('eps',[4,1]);

flag = ones(length(posax),1);

for i = 1:length(posax)
    solEps = vpasolve([dat2.ZPL(:,i) - eZPL(eps) == 0, dat2.deltaEs(:,i) - eDeltaEs(eps) == 0],eps);
    if isempty(solEps.eps1)
        solEps.eps1 = 0;
        solEps.eps2 = 0;
        solEps.eps3 = 0;
        solEps.eps4 = 0;
        flag(i) = 0;
    end
    solEpsEs(:,i) = [solEps.eps1;solEps.eps2;solEps.eps3;solEps.eps4];
    solEps = vpasolve([dat2.ZPL(:,i) - eZPL(eps) == 0, dat2.deltaGs(:,i) - eDeltaGs(eps) == 0],eps);
     if isempty(solEps.eps1)
        solEps.eps1 = 0;
        solEps.eps2 = 0;
        solEps.eps3 = 0;
        solEps.eps4 = 0;
        flag(i) = 0;
    end
    solEpsGs(:,i) = [solEps.eps1;solEps.eps2;solEps.eps3;solEps.eps4];
end

solEpsAvg = (solEpsEs + solEpsGs) / 2;
solEpsAvg(4,:) = zeros(1,length(posax));

startPoints = [solEpsAvg;zeros(2,length(posax))];

end