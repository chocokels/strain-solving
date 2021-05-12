load('peaks');
load('parameters')
load('fitPeaks.mat');
load('strainFunctions.mat');

sigma = 2;  % Peak error value
webIters = 4;  % Number of iterations to find error ranges
limNum = 2;  % Identifies which confidence interval to find (2 = 90% confidence interval)
numScans = 23;

peaks = centers * 1000;
peaksErr = zeros(4,4,26) + sigma;

confInts = [0.50 0.90 0.95 0.99];
limit = chi2inv(confInts,16);

% startPoints = getStartPoints(peaks,posax,dat0,param);  % Optional first step
startPoints = zeros(6,numScans);

[sol,chi2,solPeaks] = singleMinimizeChi(peaks,peaksErr,startPoints,dat0,param); 

ranges = zeros(6,2,numScans);
webStep = 0.5e-6;
webAx = linspace(-8e-5,8e-5,801);

for scanNum = 1:numScans
    ranges(:,:,scanNum) = web(sol(:,scanNum),peaks(:,:,scanNum),1:6,webStep,webAx,dat0,param,sigma,limit(limNum),sumSqPeaks,webIters);
end

save('solvedData.mat','sol','solPeaks','chi2','peaksErr','ranges');