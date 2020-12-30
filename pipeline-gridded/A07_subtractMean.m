close all;
clear;

%% Load data

windowSize = <PY:WINDOW_SIZE>;
minNumberOfObs = 20;

% Load mean field coefficients (learned from non-hurricane observations)
load(['./Results/meanField_',num2str(windowSize),'_',...
    num2str(minNumberOfObs),'.mat']);

% Load temperature data
% inprefix='./Data/gridTempProfHurricane_'
inprefix = '<PY:GRID_TEMP_FN>';
load([inprefix,num2str(windowSize),'_',num2str(minNumberOfObs),'.mat']);

%% Subtract mean

profLatAggrSelRounded = roundHalf(profLatAggrSel);
profLongAggrSelRounded = roundHalf(profLongAggrSel);

nProf = length(profLatAggrSelRounded);
gridTempRes = zeros(size(gridTempProf));

nPresGrid = length(presGrid);

progInterval = 20000;

parfor_progress(ceil(nProf/progInterval));

for iProf = 1:nProf
    profLat = profLatAggrSelRounded(iProf);
    profLong = profLongAggrSelRounded(iProf);
    
    iLat = find(latGrid(1,:) == profLat);
    iLong = find(longGrid(:,1) == profLong);
    
    for iPresGrid = 1:nPresGrid
        betaTemp = betaGrid(iLong,iLat,:,iPresGrid);
        
        yearDayRatio = fromJulDayToYearDay(profJulDayAggrSel(iProf))/yearLength(profJulDayAggrSel(iProf));

        tempProfHat = betaTemp(1) ...
                            + betaTemp(2) * sin(2*pi*1*yearDayRatio) + betaTemp(3) * cos(2*pi*1*yearDayRatio) ...
                            + betaTemp(4) * sin(2*pi*2*yearDayRatio) + betaTemp(5) * cos(2*pi*2*yearDayRatio) ...
                            + betaTemp(6) * sin(2*pi*3*yearDayRatio) + betaTemp(7) * cos(2*pi*3*yearDayRatio) ...
                            + betaTemp(8) * sin(2*pi*4*yearDayRatio) + betaTemp(9) * cos(2*pi*4*yearDayRatio) ...
                            + betaTemp(10) * sin(2*pi*5*yearDayRatio) + betaTemp(11) * cos(2*pi*5*yearDayRatio) ...
                            + betaTemp(12) * sin(2*pi*6*yearDayRatio) + betaTemp(13) * cos(2*pi*6*yearDayRatio) ...
                            + betaTemp(19) * (profJulDayAggrSel(iProf) - midJulDay) ...
                            + betaTemp(20) * (profJulDayAggrSel(iProf) - midJulDay)^2;
        
        gridTempRes(iProf, iPresGrid) = gridTempProf(iProf, iPresGrid) - tempProfHat;
        
    end
    if mod(iProf,progInterval) == 0
        parfor_progress;
    end
    
end

parfor_progress(0);

% outprefix = './Data/gridTempResHurricane_';
outprefix = '<PY:RES_TEMP_FN>';
save([outprefix,num2str(windowSize),'_',num2str(minNumberOfObs),'.mat'],...
    'profLatAggrSel','profLongAggrSel','profYearAggrSel',...
    'profJulDayAggrSel','profFloatIDAggrSel','profCycleNumberAggrSel',...
    'gridTempRes','intStart','intEnd','presGrid', '-v7.3');

exit;
