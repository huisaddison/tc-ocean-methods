function val = negLogLikSpaceTimeExpGeom_vec(params,profLatAggr,profLongAggr,profJulDayAggr,intTempResAggr)

thetas = exp(params(1));
thetaLat = exp(params(2));
thetaLong = exp(params(3));
thetat = exp(params(4));
sigma = exp(params(5));

disp([thetas,thetaLat,thetaLong,thetat,sigma]);

largeVal = 1e10;

if ~(all([thetas,thetaLat,thetaLong,thetat,sigma] < largeVal))
    val = NaN;
    disp(val);
    return;
end

if ~(all([thetas,thetaLat,thetaLong,thetat,sigma] > eps))
    val = NaN;
    disp(val);
    return;
end

nYear = size(intTempResAggr,2);

val = 0;

for iYear = 1:nYear
    
    profLatYear = profLatAggr{iYear};
    profLongYear = profLongAggr{iYear};
    profJulDayYear = profJulDayAggr{iYear};
    intTempResYear = intTempResAggr{iYear};
    
    nRes = length(intTempResYear);
    
    covObs = spaceTimeCovarianceExpGeom_vec(profLatYear,profLongYear,profJulDayYear,thetas,thetaLat,thetaLong,thetat);
    
%     covObs = zeros(nRes,nRes);
%     for i = 1:nRes
%         for j = 1:nRes
%             covObs(i,j) = spaceTimeCovarianceExpGeom(profLatYear(i),profLongYear(i),profJulDayYear(i),profLatYear(j),profLongYear(j),profJulDayYear(j),thetas,thetaLat,thetaLong,thetat);
%         end
%     end
    
    val = val + sum(log(eig(covObs + sigma^2*eye(nRes)))) + (intTempResYear)'*((covObs + sigma^2*eye(nRes))\intTempResYear) + log(2*pi)*nRes;
    
end

val = 0.5*val;

disp(val);