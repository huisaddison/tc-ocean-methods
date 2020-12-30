
close all;
clear;

load(['./Data/gridTempProf.mat']);

windowSize = <PY:WINDOW_SIZE>;
minNumberOfObs = 20;

load(['./Data/dataMask_',num2str(windowSize), '_',num2str(minNumberOfObs),'.mat']);

mask = dataMask;
profLatAggrRounded = roundHalf(profLatAggrSel);
profLongAggrRounded = roundHalf(profLongAggrSel);

nProf = length(profLatAggrRounded);
keep = zeros(1,nProf);

for iProf = 1:nProf
    latIdx = find(latGrid(1,:) == profLatAggrRounded(iProf));
    longIdx = find(longGrid(:,1) == profLongAggrRounded(iProf));
    
    keep(iProf) = ~isnan(mask(longIdx,latIdx));
end

disp(sum(keep));
disp(nProf);
disp(sum(keep)/nProf);

keep = logical(keep);

profLatAggrSel = profLatAggrSel(keep);
profLongAggrSel = profLongAggrSel(keep);
profYearAggrSel = profYearAggrSel(keep);
profJulDayAggrSel = profJulDayAggrSel(keep);
profFloatIDAggrSel = profFloatIDAggrSel(keep);
profCycleNumberAggrSel = profCycleNumberAggrSel(keep);
gridTempProf = gridTempProf(keep);

save(['./Data/gridTempProfFiltered_',num2str(windowSize), ...
        '_', num2str(minNumberOfObs),'.mat'],...
    'profLatAggrSel','profLongAggrSel','profYearAggrSel',...
    'profJulDayAggrSel','profFloatIDAggrSel','profCycleNumberAggrSel',...
    'gridTempProf','intStart','intEnd','-v7.3');

exit;
