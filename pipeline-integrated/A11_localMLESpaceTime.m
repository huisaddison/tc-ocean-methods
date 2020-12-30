close all;
clear;

windowSize = <PY:WINDOW_SIZE>;
minNumberOfObs = 20;
% region = '_WesternPacific';
region = '<PY:OCEAN_BASIN>';

month = <PY:CENTER_MONTH>;

startYear = <PY:START_YEAR>;
endYear = <PY:END_YEAR>;

[latGrid,longGrid] = <PY:OB_MESHGRID>;

nGrid = numel(latGrid);
nYear = endYear - startYear + 1;

thetasOpt = zeros(size(latGrid));
thetaLatOpt = zeros(size(latGrid));
thetaLongOpt = zeros(size(latGrid));
thetatOpt = zeros(size(latGrid));
sigmaOpt = zeros(size(latGrid));
nll = zeros(size(latGrid));
nResGrid = zeros(size(latGrid));

% Amended to smaller window size
windowSizeGP = 5;

% Discard previous iterIdx, if it exists
fileName = ['iterIdxLocalMLESpaceTime'...
    '_',num2str(windowSize),'_',num2str(minNumberOfObs),'_',...
    num2str(windowSizeGP),'_',num2str(month,'%02d'),'_',...
    num2str(startYear),'_',num2str(endYear),region,'.txt'];
fileID = fopen(fileName,'w');
fclose(fileID);

parpool(<PY:N_PARPOOL>)

tic;
parfor iGrid = 1:nGrid
    
    fileID = fopen(fileName,'a');
    fprintf(fileID,'%d \n',iGrid);
    fclose(fileID);
    
    predLat = latGrid(iGrid);
    predLong = longGrid(iGrid);

    latMin = predLat - windowSizeGP;
    latMax = predLat + windowSizeGP;
    longMin = predLong - windowSizeGP;
    longMax = predLong + windowSizeGP;

    profLatAggr = cell(1,nYear);
    profLongAggr = cell(1,nYear);
    profJulDayAggr = cell(1,nYear);
    tempResAggr = cell(1,nYear);

    for iYear = startYear:endYear

        if strcmp(region, '_AllBasins')
            S = load(['./Data/Extended/gridTempRes_',num2str(windowSize),...
                '_',num2str(minNumberOfObs),'_Month_',num2str(month,'%02d'),...
                '_',num2str(iYear), '_extended.mat']);
        else 
            S = load(['./Data/Extended/gridTempRes_',num2str(windowSize),...
                '_',num2str(minNumberOfObs),'_Month_',num2str(month,'%02d'),...
                '_',num2str(iYear),'_extended_filtered',region,'.mat']);
        end

        
        profLat3Months = S.profLatAggr3Months;
        profLong3Months = S.profLongAggr3Months;
        profJulDay3Months = S.profJulDayAggr3Months;
        tempRes3Months = S.gridTempRes3Months(:)';
        
        % Enable wrap around by duplicating boundary data
        leftBoundaryIdx = find(profLong3Months <= 20 + windowSizeGP);
        rightBoundaryIdx = find(profLong3Months >= 380 - windowSizeGP);
        profLong3Months = [profLong3Months ...
            profLong3Months(leftBoundaryIdx) + 360 ...
            profLong3Months(rightBoundaryIdx) - 360];
        profLat3Months = [profLat3Months ...
            profLat3Months(leftBoundaryIdx) ...
            profLat3Months(rightBoundaryIdx)];
        profJulDay3Months = [profJulDay3Months ...
            profJulDay3Months(leftBoundaryIdx) ...
            profJulDay3Months(rightBoundaryIdx)];
        tempRes3Months = [tempRes3Months ...
            tempRes3Months(leftBoundaryIdx) ...
            tempRes3Months(rightBoundaryIdx)];

        idx = find(profLat3Months > latMin ...
            & profLat3Months < latMax ...
            & profLong3Months > longMin ...
            & profLong3Months < longMax);

        profLatAggr{iYear-startYear+1} = profLat3Months(idx)';
        profLongAggr{iYear-startYear+1} = profLong3Months(idx)';
        profJulDayAggr{iYear-startYear+1} = profJulDay3Months(idx)';
        tempResAggr{iYear-startYear+1} = tempRes3Months(idx)';

    end
    
    nResGrid(iGrid) = sum(cellfun(@length,tempResAggr));
    
    if nResGrid(iGrid) == 0 % No observations in the window
        
        thetasOpt(iGrid) = NaN;
        thetaLatOpt(iGrid) = NaN;
        thetaLongOpt(iGrid) = NaN;
        thetatOpt(iGrid) = NaN;
        sigmaOpt(iGrid) = NaN;
        nll(iGrid) = NaN;
        
        continue;
    end

    try
        fun = @(params) negLogLikSpaceTimeExpGeom_vec(params,...
            profLatAggr,profLongAggr,profJulDayAggr,tempResAggr);
        
        logThetasInit = log(400^2);
        logThetaLatInit = log(5);
        logThetaLongInit = log(5);
        logThetatInit = log(5);
        logSigmaInit = log(100);
        
        opts = optimoptions(@fminunc,'Algorithm','quasi-newton',...
            'MaxFunctionEvaluations',1000);

        [paramOpt,nll(iGrid)] = fminunc(fun,[logThetasInit, logThetaLatInit,...
            logThetaLongInit, logThetatInit, logSigmaInit],opts);
        
        thetasOpt(iGrid) = exp(paramOpt(1));
        thetaLatOpt(iGrid) = exp(paramOpt(2));
        thetaLongOpt(iGrid) = exp(paramOpt(3));
        thetatOpt(iGrid) = exp(paramOpt(4));
        sigmaOpt(iGrid) = exp(paramOpt(5));

    catch
        warning('Optimization failed!');

        thetasOpt(iGrid) = NaN;
        thetaLatOpt(iGrid) = NaN;
        thetaLongOpt(iGrid) = NaN;
        thetatOpt(iGrid) = NaN;
        sigmaOpt(iGrid) = NaN;
        nll(iGrid) = NaN;
    end
    
end
toc;

save(['./Results/localMLESpaceTime_Depth_',...
        num2str(windowSize),'_',num2str(minNumberOfObs),'_',...
        num2str(windowSizeGP),'_',num2str(month,'%02d'),'_',...
        num2str(startYear),'_',num2str(endYear),region,'.mat'],...
    'latGrid','longGrid','thetasOpt','thetaLatOpt','thetaLongOpt',...
    'thetatOpt','sigmaOpt','nll','nResGrid','-v7.3');
exit;
