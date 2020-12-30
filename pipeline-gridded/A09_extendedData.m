close all;
clear;

windowSize = <PY:WINDOW_SIZE>;
minNumberOfObs = 20;

load_prefix = ['./Data/Monthly/gridTempResNonHurricane_',...
    num2str(windowSize),'_',num2str(minNumberOfObs),...
    '_Month_'];

for iYear = <PY:START_YEAR>:<PY:END_YEAR>
    for iMonth = 1:12
        % Data for iMonth-1
        if iMonth == 1
            if iYear == <PY:START_YEAR>
                S1.profLatAggrMonth = [];
                S1.profLongAggrMonth = [];
                S1.profJulDayAggrMonth = [];
                S1.profFloatIDAggrMonth = [];
                S1.profCycleNumberAggrMonth = [];
                S1.gridTempResMonth = [];
            else
                S1 = load([load_prefix,num2str(12,'%02d'),...
                    '_',num2str(iYear-1),'.mat']);
            end
        else
            S1 = load([load_prefix,num2str(iMonth-1,'%02d'),...
                '_',num2str(iYear),'.mat']);
        end
        
        % Data for iMonth        
        S2 = load([load_prefix,num2str(iMonth,'%02d'),...
            '_',num2str(iYear),'.mat']);
        
        % Data for iMonth+1
        if iMonth == 12
            if iYear == <PY:END_YEAR>
                S3.profLatAggrMonth = [];
                S3.profLongAggrMonth = [];
                S3.profJulDayAggrMonth = [];
                S3.profFloatIDAggrMonth = [];
                S3.gridTempResMonth = [];
            else
                S3 = load([load_prefix,num2str(1,'%02d'),...
                    '_',num2str(iYear+1),'.mat']);
            end
        else
            S3 = load([load_prefix,num2str(iMonth+1,'%02d'),...
                '_',num2str(iYear),'.mat']);
        end
        
        profLatAggr3Months = [S1.profLatAggrMonth S2.profLatAggrMonth...
            S3.profLatAggrMonth];
        profLongAggr3Months = [S1.profLongAggrMonth S2.profLongAggrMonth...
            S3.profLongAggrMonth];
        profJulDayAggr3Months =  [S1.profJulDayAggrMonth ...
            S2.profJulDayAggrMonth S3.profJulDayAggrMonth];
        profFloatIDAggr3Months = [S1.profFloatIDAggrMonth ...
            S2.profFloatIDAggrMonth S3.profFloatIDAggrMonth];
        profCycleNumberAggr3Months = [S1.profCycleNumberAggrMonth ...
            S2.profCycleNumberAggrMonth S3.profCycleNumberAggrMonth];
        gridTempRes3Months = [S1.gridTempResMonth; ...
            S2.gridTempResMonth; S3.gridTempResMonth];
        intStart = S2.intStart;
        intEnd   = S2.intEnd;
        presGrid = S2.presGrid;
        
        save(['./Data/Extended/gridTempRes_',num2str(windowSize),...
            '_',num2str(minNumberOfObs),'_Month_',...
            num2str(iMonth,'%02d'),'_',num2str(iYear),'_extended.mat'],...
        'gridTempRes3Months','profLatAggr3Months',...
        'profLongAggr3Months','profJulDayAggr3Months',... 
        'profFloatIDAggr3Months','profCycleNumberAggr3Months',...
        'intStart','intEnd','presGrid');
        
    end
end
exit;
