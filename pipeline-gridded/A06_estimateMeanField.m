close all;
clear;

windowSize = <PY:WINDOW_SIZE>;

cp0 = 3989.244;
rho0 = 1030;

formatIn = 'dd-mmm-yyyy';
startDay = datenum('01-Jan-<PY:START_YEAR>', formatIn);
endDay = datenum('01-Jan-<PY:END_YEAR+1>', formatIn);
minNumberOfObs = 20;
midJulDay = (endDay - startDay) / 2 + startDay


%% Load integrated profiles

load(['./Data/gridTempProfNonHurricane_',num2str(windowSize),'_',num2str(minNumberOfObs),'.mat']);
nPresGrid = length(presGrid);

%% Enable wrap around by duplicating boundary data

leftBoundaryIdx = find(profLongAggrSel <= 20 + windowSize);
rightBoundaryIdx = find(profLongAggrSel >= 380 - windowSize);
profLongAggrSel = [profLongAggrSel profLongAggrSel(leftBoundaryIdx) + 360 profLongAggrSel(rightBoundaryIdx) - 360];
profLatAggrSel = [profLatAggrSel profLatAggrSel(leftBoundaryIdx) profLatAggrSel(rightBoundaryIdx)];
profJulDayAggrSel = [profJulDayAggrSel profJulDayAggrSel(leftBoundaryIdx) profJulDayAggrSel(rightBoundaryIdx)];
gridTempProf = [gridTempProf; gridTempProf(leftBoundaryIdx, :); gridTempProf(rightBoundaryIdx, :)];

%% Calculate mean field using a moving window

[latGrid,longGrid] = meshgrid(linspace(-89.5,89.5,180),linspace(20.5,379.5,360));

nGrid = numel(latGrid);

%% Modify to do all at once
betaGrid = zeros([size(latGrid),20,nPresGrid]);

load(['./Data/dataMask_',num2str(windowSize),'_',num2str(minNumberOfObs),'.mat']);

mask = dataMask;

tic;

parfor_progress(nGrid);

for iGrid = 1:nGrid
    latSel = latGrid(iGrid);
    longSel = longGrid(iGrid);

    latMin = latSel - windowSize;
    latMax = latSel + windowSize;
    longMin = longSel - windowSize;
    longMax = longSel + windowSize;
    
    idx = find(profLatAggrSel > latMin & profLatAggrSel < latMax & profLongAggrSel > longMin & profLongAggrSel < longMax);
    
    latIdx = find(latGrid(1,:) == latSel);
    longIdx = find(longGrid(:,1) == longSel);    
    goodPixel = ~isnan(mask(longIdx,latIdx));
    
    % Need at least 20 data points to estimate the regression coefficients, also do not compute mean if outside land/datamask
    if length(idx) < 20 || ~goodPixel
        parfor_progress();
        continue; 
    end
    
    profJulDayAggrWindow = profJulDayAggrSel(idx)';
    profLatAggrWindow = profLatAggrSel(idx)';
    profLongAggrWindow = profLongAggrSel(idx)';
    profYearDayAggrWindow = fromJulDayToYearDay(profJulDayAggrWindow);
    profYearLengthAggrWindow = yearLength(profJulDayAggrWindow);

    % Iterate and fit over all depths
    for iPresGrid = 1:nPresGrid
        gridTempProfWindow = gridTempProf(idx, iPresGrid);
        
        %% NoTrend
        XWindow = [ones(length(profJulDayAggrWindow),1) ...
           sin(2*pi*1*profYearDayAggrWindow./profYearLengthAggrWindow) cos(2*pi*1*profYearDayAggrWindow./profYearLengthAggrWindow) ...
           sin(2*pi*2*profYearDayAggrWindow./profYearLengthAggrWindow) cos(2*pi*2*profYearDayAggrWindow./profYearLengthAggrWindow) ...
           sin(2*pi*3*profYearDayAggrWindow./profYearLengthAggrWindow) cos(2*pi*3*profYearDayAggrWindow./profYearLengthAggrWindow) ...
           sin(2*pi*4*profYearDayAggrWindow./profYearLengthAggrWindow) cos(2*pi*4*profYearDayAggrWindow./profYearLengthAggrWindow) ...
           sin(2*pi*5*profYearDayAggrWindow./profYearLengthAggrWindow) cos(2*pi*5*profYearDayAggrWindow./profYearLengthAggrWindow) ...
           sin(2*pi*6*profYearDayAggrWindow./profYearLengthAggrWindow) cos(2*pi*6*profYearDayAggrWindow./profYearLengthAggrWindow) ...
           (profLatAggrWindow-latSel) (profLongAggrWindow-longSel) (profLatAggrWindow-latSel).*(profLongAggrWindow-longSel) ...
           (profLatAggrWindow-latSel).^2 (profLongAggrWindow-longSel).^2];
            
        betaWindow = XWindow\gridTempProfWindow;
        betaWindow = [betaWindow; 0; 0]; % Put zero coefficients in place of the time trend terms
        [iGridSub1,iGridSub2] = ind2sub(size(latGrid),iGrid);
        betaGrid(iGridSub1,iGridSub2,:,iPresGrid) = betaWindow;
    end
    parfor_progress;
    
end

parfor_progress(0);

toc;

save(['./Results/meanField_',num2str(windowSize),'_',num2str(minNumberOfObs),'.mat'],'betaGrid','latGrid','longGrid','midJulDay','presGrid', '-v7.3');

exit;
