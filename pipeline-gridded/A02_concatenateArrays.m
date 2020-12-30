prefix = './Data/gridTempProf';

A = load([prefix, '_2007_2010.mat']);
B = load([prefix, '_2011_2014.mat']);
C = load([prefix, '_2015_2016.mat']);
D = load([prefix, '_2017_2018.mat']);
gridTempProf = cat(1, A.gridTempProf, B.gridTempProf, C.gridTempProf, D.gridTempProf);
profLatAggrSel = cat(2, A.profLatAggrSel, B.profLatAggrSel, C.profLatAggrSel, D.profLatAggrSel);
profLongAggrSel = cat(2, A.profLongAggrSel, B.profLongAggrSel, C.profLongAggrSel, D.profLongAggrSel);
profYearAggrSel = cat(2, A.profYearAggrSel, B.profYearAggrSel, C.profYearAggrSel, D.profYearAggrSel);
profJulDayAggrSel = cat(2, A.profJulDayAggrSel, B.profJulDayAggrSel, C.profJulDayAggrSel, D.profJulDayAggrSel);
profFloatIDAggrSel = cat(2, A.profFloatIDAggrSel, B.profFloatIDAggrSel, C.profFloatIDAggrSel, D.profFloatIDAggrSel);
profCycleNumberAggrSel = cat(2, A.profCycleNumberAggrSel, B.profCycleNumberAggrSel, C.profCycleNumberAggrSel, D.profCycleNumberAggrSel);
presGrid=A.presGrid;
intStart=A.intStart;
intEnd=A.intEnd;
save(['./Data/gridTempProf.mat'], 'profLatAggrSel','profLongAggrSel','profYearAggrSel','profJulDayAggrSel','profFloatIDAggrSel','profCycleNumberAggrSel','gridTempProf','intStart','intEnd','presGrid', '-v7.3');

exit;
