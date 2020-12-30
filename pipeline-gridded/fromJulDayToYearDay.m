function r = fromJulDayToYearDay(julDay)    

tempDateVec = datevec(julDay);
r = datenum(tempDateVec) - datenum([tempDateVec(:,1) repmat([1 1 0 0 0],length(julDay),1)]);