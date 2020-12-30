load('OBM.mat')

NotPacific = 1 - pacificmask;
NorthAtlantic = zeros(180, 360);
NorthAtlantic(90:160, 260:360) = 1;

[earthLatT, earthLongT] = meshgrid(1:180, 1:360);
earthLat = flipud(earthLatT');
earthLong = earthLongT';
