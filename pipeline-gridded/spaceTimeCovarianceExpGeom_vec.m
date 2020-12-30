function r = spaceTimeCovarianceExpGeom_vec(lat,long,time,thetas,thetaLat,thetaLong,thetat)

distSpaceTime = sqrt(((lat - lat')./thetaLat).^2 + ((long - long')./thetaLong).^2 + ((time - time')./thetat).^2);

r = thetas .* exp(-distSpaceTime);