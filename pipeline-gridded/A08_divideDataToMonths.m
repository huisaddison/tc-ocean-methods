close all;
clear;

windowSize = <PY:WINDOW_SIZE>;
minNumberOfObs = 20;

load(['./Data/gridTempResNonHurricane_',num2str(windowSize),'_',...
    num2str(minNumberOfObs),'.mat']);


nProf = length(profLatAggrSel);

profMonthAggrSel = zeros(1,nProf);
for iProf = 1:nProf
    tempv = datevec(profJulDayAggrSel(iProf));
    profMonthAggrSel(iProf) = tempv(2);
end

for iYear = <PY:START_YEAR>:<PY:END_YEAR>
    for iMonth = 1:12
        idx = (profYearAggrSel == iYear & profMonthAggrSel == iMonth);
        
        gridTempResMonth = gridTempRes(idx, :);
        profLatAggrMonth = profLatAggrSel(idx);
        profLongAggrMonth = profLongAggrSel(idx);
        profFloatIDAggrMonth = profFloatIDAggrSel(idx);
        profJulDayAggrMonth = profJulDayAggrSel(idx);
        profCycleNumberAggrMonth = profCycleNumberAggrSel(idx);
        
        save(['./Data/Monthly/gridTempResNonHurricane_',...
            num2str(windowSize),'_',num2str(minNumberOfObs),...
            '_Month_',num2str(iMonth,'%02d'),'_',num2str(iYear),'.mat'],...
        'gridTempResMonth','profLatAggrMonth','profLongAggrMonth',...
        'profFloatIDAggrMonth','profJulDayAggrMonth',...
        'profCycleNumberAggrMonth','intStart','intEnd','presGrid');
        
        disp(sum(idx));
    end
end
exit;
