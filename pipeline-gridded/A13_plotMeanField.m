close all;
clear;


%% Plotting

cp0 = 3989.244;
rho0 = 1030;

windowSize = <PY:WINDOW_SIZE>;
FigFolder = 'Figures_Debug';
load(['./Results/meanField_',num2str(windowSize),'_20.mat']);
load(['./Data/dataMask_',num2str(windowSize),'_20.mat']);

mask = dataMask;

for idx = 1:length(presGrid)
    depth = presGrid(idx);
    disp(depth);
    %% Plot OHC trend
    figure;
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off');
    unitsScaling = cp0 * rho0 / (24*60*60);
    surfm(latGrid,longGrid,unitsScaling*betaGrid(:,:,19,idx));
    load coast;
    plotm(lat,long,'k');
    h = colorbar;
    caxis([-4,4]);
    cLims = caxis;
    colormap(darkb2r(cLims(1),cLims(2)));
    drawnow;
    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 22.5 15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    print('-depsc2',['./',FigFolder,'/OHC_trend_',...
        num2str(depth, '%03d'),'.eps']);

    %% Plot OHC trend 2

    figure;
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off');
    surfm(latGrid,longGrid,betaGrid(:,:,20,idx));
    load coast;
    plotm(lat,long,'k');
    h = colorbar;
    temp = betaGrid(:,:,20,idx);
    caxis([quantile(temp(:),0.01),-quantile(temp(:),0.01)]);
    cLims = caxis;
    colormap(darkb2r(cLims(1),cLims(2)));
    drawnow;
    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 22.5 15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    print('-depsc2',['./',FigFolder,'/OHC_trend2_',...
        num2str(depth, '%03d'),'.eps']);

    %% Plot beta0
    figure;
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off');
    surfm(latGrid,longGrid,mask.*betaGrid(:,:,1,idx));
    load coast;
    plotm(lat,long,'k');
    h = colorbar;
    colormap(jet(100));
    colorbar;
    caxis([-1e3,0.5e4]);

    drawnow;

    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 22.5 15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    print('-depsc2',['./',FigFolder,'/beta0_',...
        num2str(depth, '%03d'),'.eps']);

    %% Plot sin harmonics

    for iHar = 1:6

        figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');
        sinHar = iHar;
        surfm(latGrid,longGrid,betaGrid(:,:,2*sinHar,idx));
        load coast;
        plotm(lat,long,'k');
        h = colorbar;
        temp = betaGrid(:,:,2*sinHar,idx);
        caxis([quantile(temp(:),0.01),-quantile(temp(:),0.01)]);
        cLims = caxis;
        colormap(darkb2r(cLims(1),cLims(2)));
        drawnow;
        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 22.5 15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        print('-depsc2',['./',FigFolder,'/betaSin_',...
            num2str(sinHar), '_',...
            num2str(depth, '%03d'),'.eps']);

    end

    %% Plot cos harmonics

    for iHar = 1:6
        figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');
        cosHar = iHar;
        surfm(latGrid,longGrid,betaGrid(:,:,2*cosHar+1,idx));
        load coast;
        plotm(lat,long,'k');
        h = colorbar;
        temp = betaGrid(:,:,2*cosHar,idx);
        caxis([quantile(temp(:),0.01),-quantile(temp(:),0.01)]);
        cLims = caxis;
        colormap(darkb2r(cLims(1),cLims(2)));
        drawnow;
        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 22.5 15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        print('-depsc2',['./',FigFolder,'/betaCos_',...
            num2str(cosHar), '_',...
            num2str(depth, '%03d'),'.eps']);
    end

    %% Plot leading harmonic amplitude
    figure;
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off');
    surfm(latGrid,longGrid,mask.*sqrt(betaGrid(:,:,2,idx).^2 ...
        + betaGrid(:,:,3,idx).^2));

    load coast;
    plotm(lat,long,'k');
    colormap(jet(100));
    colorbar;
    caxis([0,100]);
    drawnow;
    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 22.5 15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    print('-depsc2',['./',FigFolder,'/amplitude_',...
        num2str(depth, '%03d'),'.eps']);

    %% Plot leading harmonic phase shift
    figure;
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off');
    surfm(latGrid,longGrid,mask.*atan2(betaGrid(:,:,3,idx),...
        betaGrid(:,:,2,idx)));
    load coast;
    plotm(lat,long,'k');
    colormap([parula(100); flipud(parula(100))]);
    colorbar;
    caxis([-pi,pi]);
    drawnow;
    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 22.5 15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    print('-depsc2',['./',FigFolder,'/phaseShift_',...
        num2str(depth, '%03d'),'.eps']);

    %% Plot betaLat
    figure;
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off');
    surfm(latGrid,longGrid,betaGrid(:,:,14,idx));
    load coast;
    plotm(lat,long,'k');
    h = colorbar;
    caxis([-40,40]);
    cLims = caxis;
    colormap(darkb2r(cLims(1),cLims(2)));
    drawnow;
    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 22.5 15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    print('-depsc2',['./',FigFolder,'/betaLat_',...
        num2str(depth, '%03d'),'.eps']);

    %% Plot betaLong
    figure;
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off');
    surfm(latGrid,longGrid,betaGrid(:,:,15,idx));
    load coast;
    plotm(lat,long,'k');
    h = colorbar;
    caxis([-20,20]);
    cLims = caxis;
    colormap(darkb2r(cLims(1),cLims(2)));
    drawnow;
    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 22.5 15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    print('-depsc2',['./',FigFolder,'/betaLong_',...
        num2str(depth, '%03d'),'.eps']);

    %% Plot betaLatLong

    figure;
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off'); 
    surfm(latGrid,longGrid,betaGrid(:,:,16,idx)); 
    load coast;
    plotm(lat,long,'k'); 
    h = colorbar; 
    caxis([-5,5]); 
    cLims = caxis;
    colormap(darkb2r(cLims(1),cLims(2)));
    drawnow;
    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 22.5 15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    print('-depsc2',['./',FigFolder,'/betaLatLong_',...
        num2str(depth, '%03d'),'.eps']);

    %% Plot betaLat2 
    figure;
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off'); 
    surfm(latGrid,longGrid,betaGrid(:,:,17,idx)); 
    load coast;
    plotm(lat,long,'k'); 
    h = colorbar; 
    caxis([-5,5]); 
    cLims = caxis;
    colormap(darkb2r(cLims(1),cLims(2))); 
    drawnow; 
    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 22.5 15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    print('-depsc2',['./',FigFolder,'/betaLat2_',...
        num2str(depth, '%03d'),'.eps']);

    %% Plot betaLong2 
    figure;
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off'); 
    surfm(latGrid,longGrid,betaGrid(:,:,18,idx)); 
    load coast;
    plotm(lat,long,'k'); 
    h = colorbar; 
    caxis([-5,5]); 
    cLims = caxis;
    colormap(darkb2r(cLims(1),cLims(2))); 
    drawnow; 
    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 22.5 15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    print('-depsc2',['./',FigFolder,'/betaLong2_',...
        num2str(depth, '%03d'),'.eps']);

end
