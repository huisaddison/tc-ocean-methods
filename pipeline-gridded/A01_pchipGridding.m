close all;
clear;

PATH_TO_GSW='/home/ahu2/GSW/';
PATH_TO_AGGR_DATA='<PY:DATA_LOC>'
addpath(genpath(PATH_TO_GSW));


%% Load data
years = '<PY:YEARS>';
load([PATH_TO_AGGR_DATA, 'Argo_data_aggr_', years, '.mat']);

%% Filter out duplicate profiles
[C,ia,ic] = unique([profLatAggr',profLongAggr',profJulDayAggr'],'rows');

profPresAggrGood = profPresAggr(ia);
profTempAggrGood = profTempAggr(ia);
profPsalAggrGood = profPsalAggr(ia);
profLatAggrGood = profLatAggr(ia);
profLongAggrGood = profLongAggr(ia);
profYearAggrGood = profYearAggr(ia);
profJulDayAggrGood = profJulDayAggr(ia);
profFloatIDAggrGood = profFloatIDAggr(ia);
profCycleNumberAggrGood = profCycleNumberAggr(ia);

%% Starting pressure and end pressure histograms

startPres = cellfun(@min,profPresAggrGood);
endPres = cellfun(@max,profPresAggrGood);

%% Profile selection based on start and end pressure
intStart = <PY:GRID_LOWER>;
intEnd   = <PY:GRID_UPPER>;
selIdx = (startPres >=0 & startPres <= intStart & endPres >= intEnd);

profPresAggrSel = profPresAggrGood(selIdx);
profTempAggrSel = profTempAggrGood(selIdx);
profPsalAggrSel = profPsalAggrGood(selIdx);
profLatAggrSel = profLatAggrGood(selIdx);
profLongAggrSel = profLongAggrGood(selIdx);
profYearAggrSel = profYearAggrGood(selIdx);
profJulDayAggrSel = profJulDayAggrGood(selIdx);
profFloatIDAggrSel = profFloatIDAggrGood(selIdx);
profCycleNumberAggrSel = profCycleNumberAggrGood(selIdx);


%% Compute absolute salinity and conservative and potential temperature, you'll need to have the GSW toolbox in Matlab path to run this section, see http://www.teos-10.org/software.htm

% Convert longitude from 20-380 range to 0-360 range
profLongAggrSelTemp = (profLongAggrSel > 360).*(profLongAggrSel - 360) + (profLongAggrSel <= 360).*profLongAggrSel;

profAbsSalAggrSel = cellfun(@gsw_SA_from_SP,profPsalAggrSel, profPresAggrSel,...
    num2cell(profLongAggrSelTemp), num2cell(profLatAggrSel), 'UniformOutput', 0);
profPotTempAggrSel = cellfun(@gsw_pt_from_t,profAbsSalAggrSel,profTempAggrSel,...
    profPresAggrSel,'UniformOutput',0);


presGrid = <PY:GRID_LOWER>:<PY:GRID_STRIDE>:<PY:GRID_UPPER>;
nPresGrid = length(presGrid);
nProf = length(profPresAggrSel);
gridTempProf = zeros(nProf, nPresGrid);

% parpool(<PY:N_PARPOOL>)
% parfor_progress(nProf);

for i = 1:nProf
    press=profPresAggrSel{i};
    pottemp=profPotTempAggrSel{i};
    gridTempProf(i, :)=pchip(press, pottemp, presGrid);
    % parfor_progress;
end
% parfor_progress(0);

save(['./Data/gridTempProf_',years,'.mat'],'profLatAggrSel','profLongAggrSel','profYearAggrSel','profJulDayAggrSel','profFloatIDAggrSel','profCycleNumberAggrSel','gridTempProf','intStart','intEnd','presGrid','-v7.3');

exit;
