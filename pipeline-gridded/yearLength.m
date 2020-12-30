function r = yearLength(julDay)

tempDateVec = datevec(julDay);
r = datenum(tempDateVec(:,1),12,31,24,0,0)-datenum(tempDateVec(:,1),1,1,0,0,0);