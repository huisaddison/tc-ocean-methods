close all;
clear;

%% Load data

verticalSelectionImpacted = 'Hurricane';
verticalSelectionBackground = 'NoHurs';
hprefix='<PY:GRID_TEMP_FN>';
windowSize=<PY:WINDOW_SIZE>;
windowSizeGP=5;
minNumberOfObs = 20;
month = <PY:CENTER_MONTH>;
startYear =<PY:START_YEAR>;
endYear =<PY:END_YEAR>;
region = '<PY:OCEAN_BASIN>';

if strcmp(region, '_WesternPacific')
        % WesternPacific _only_ grid 
            %[latGrid,longGrid] = meshgrid(0.5:65.5,105.5:187.5);
        latMin = 0;
        latMax = 70;
        longMin = 105;
        longMax = 188;
elseif strcmp(region, '_NorthAtlantic')
        % NorthPacific _only_ grid 
            %[latGrid,longGrid] = meshgrid(0.5:70.5,261.5:360.5);
        latMin = 0;
        latMax = 71;
        longMin = 262;
        longMax = 361;
elseif strcmp(region, '_AllBasins')
        latMin = -90;
        latMax = 90;
        longMin = 20;
        longMax = 380;
else 
        disp('Region unrecognized');
        exit;
end 

% We are estimating the background process for hurricane-impacted
% float pairs, so we load:
% - hurricane-impacted float integrated heat content
% - coefficients learned from non-hurricane-impacted floats
H = load([hprefix,num2str(windowSize), '_', num2str(minNumberOfObs),'.mat']);
B = load(['./Results/localMLESpaceTime_Depth_',...
        num2str(windowSize),'_',num2str(minNumberOfObs),'_',...
        num2str(windowSizeGP),'_',num2str(month,'%02d'),'_',...
        num2str(startYear),'_',num2str(endYear),region,'.mat']);

% Filter hurricane observations by ocean region
% Necessary as coefficients are stored in separate files by
% ocean region.

inRegion =    (H.profLatAggrSel > latMin) ...
           .* (H.profLatAggrSel < latMax) ...
           .* (H.profLongAggrSel > longMin) ...
           .* (H.profLongAggrSel < longMax) > 0;
CycleNumberReg =    H.profCycleNumberAggrSel(inRegion);
FloatIDReg =        H.profFloatIDAggrSel(inRegion);
JulDayReg =         H.profJulDayAggrSel(inRegion);
LatReg =            H.profLatAggrSel(inRegion);
LongReg =           H.profLongAggrSel(inRegion);
YearReg =           H.profYearAggrSel(inRegion);

LatRegRounded = roundHalf(LatReg);
LongRegRounded = roundHalf(LongReg);
nProf = length(CycleNumberReg);

%% Create new coefficient vectors
sigma =             zeros(1, nProf);
thetaLat =          zeros(1, nProf);
thetaLong =         zeros(1, nProf);
thetas =            zeros(1, nProf);
thetat =            zeros(1, nProf);
nll =               zeros(1, nProf);
nResGrid =          zeros(1, nProf);
interpLat =         zeros(1, nProf);
interpLong =        zeros(1, nProf);

progInterval = 2000;

parfor_progress(ceil(nProf/progInterval));


for iProf = 1:nProf
    % Iterate and store rounded
    profLat = LatRegRounded(iProf);
    profLong = LongRegRounded(iProf);
	  interpLat(iProf) = profLat;
    interpLong(iProf) = profLong;

    % Locate coefficient indices
    iLat = find(B.latGrid(1,:) == profLat);
    iLong = find(B.longGrid(:,1) == profLong);
    
    % Read off coefficients
    sigma(iProf)        = B.sigmaOpt(iLong, iLat);
    thetaLat(iProf)     = B.thetaLatOpt(iLong, iLat);
    thetaLong(iProf)    = B.thetaLongOpt(iLong, iLat);
    thetas(iProf)       = B.thetasOpt(iLong, iLat);
    thetat(iProf)       = B.thetatOpt(iLong, iLat);
    nll(iProf)          = B.nll(iLong, iLat);
    nResGrid(iProf)     = B.nResGrid(iLong, iLat);

    if mod(iProf,progInterval) == 0
        parfor_progress;
    end

end

parfor_progress(0);


%% Save results
save(['./Results/localMLESpaceTimeCoefs', meanTag,tag,verticalSelectionImpacted,'_',num2str(windowSize), '_',num2str(minNumberOfObs),'_',num2str(windowSizeGP), '_',num2str(month,'%02d'),'_',num2str(startYear),'_', num2str(endYear),region,'.mat'], 'FloatIDReg', 'CycleNumberReg', 'interpLat', 'interpLong', 'sigma', 'thetaLat', 'thetaLong', 'thetas', 'thetat', 'nll', 'nResGrid', 'interpLat', 'interpLong', '-v7.3');
