close all;
clear;

load(['./Data/gridTempProf.mat']);

[latGrid,longGrid] = meshgrid(linspace(-89.5,89.5,180),linspace(20.5,379.5,360));

%% Count number of data points for each month and year within each window

windowSize = <PY:WINDOW_SIZE>;

% Enable wrap around by duplicating boundary data
leftBoundaryIdx = find(profLongAggrSel <= 20 + windowSize);
rightBoundaryIdx = find(profLongAggrSel >= 380 - windowSize);
profLongAggrSel = [profLongAggrSel profLongAggrSel(leftBoundaryIdx) + 360 profLongAggrSel(rightBoundaryIdx) - 360];
profLatAggrSel = [profLatAggrSel profLatAggrSel(leftBoundaryIdx) profLatAggrSel(rightBoundaryIdx)];
profJulDayAggrSel = [profJulDayAggrSel profJulDayAggrSel(leftBoundaryIdx) profJulDayAggrSel(rightBoundaryIdx)];

nGrid = numel(latGrid);

yearMonthCounts = zeros([size(latGrid),10,12]);

parfor_progress(nGrid);

for iGrid = 1:nGrid
    latSel = latGrid(iGrid);
    longSel = longGrid(iGrid);

    latMin = latSel - windowSize;
    latMax = latSel + windowSize;
    longMin = longSel - windowSize;
    longMax = longSel + windowSize;
    
    idx = find(profLatAggrSel > latMin & profLatAggrSel < latMax & profLongAggrSel > longMin & profLongAggrSel < longMax);
    
    profJulDayAggrWindow = profJulDayAggrSel(idx)';    
    
    [iGridSub1,iGridSub2] = ind2sub(size(latGrid),iGrid);
    
    nProf = length(idx);
    profYearAggrWindow = zeros(1,nProf);
    profMonthAggrWindow = zeros(1,nProf);
    for iProf = 1:nProf
        temp = datevec(profJulDayAggrWindow(iProf));
        profYearAggrWindow(iProf) = temp(1);
        profMonthAggrWindow(iProf) = temp(2);
    end
    
    years = 2007:2018;
    yearMonthCountsWindow = zeros(10,12);
    for iYear=1:10
        yearMonthCountsWindow(iYear,:) = histcounts(profMonthAggrWindow(profYearAggrWindow == years(iYear)),0.5:12.5);
    end
    
    yearMonthCounts(iGridSub1,iGridSub2,:,:) = yearMonthCountsWindow;
    parfor_progress;
    
end

parfor_progress(0);

%% Form data-driven landmask

minNumberOfObs = 20;

dataMask = zeros(size(latGrid));

for iGrid = 1:nGrid
    
    [iGridSub1,iGridSub2] = ind2sub(size(latGrid),iGrid);
    
    predLat = latGrid(iGrid);
    predLong = longGrid(iGrid);
    
    yearMonthCountsWindow = squeeze(yearMonthCounts(iGridSub1,iGridSub2,:,:));
    
    dataMask(iGrid) = all(sum(yearMonthCountsWindow,1) >= minNumberOfObs);
    
end

dataMask(dataMask == 0) = NaN;

save(['./Data/dataMask_',num2str(windowSize), '_', ...
        num2str(minNumberOfObs),'.mat'],...
    'dataMask','latGrid','longGrid','-v7.3');

exit;
