load('sweepMay24N500.mat')

atindex = 501;
bs(atindex)

atrate_mf = zeros(length(wrange),length(kaprange));
atrate_sm = zeros(length(wrange),length(kaprange));
betasteps = 1000;
for wi=1:length(wrange)
    w=wrange(wi)
    for ki=1:length(kaprange)
        kap=kaprange(ki);
        rho_est  = ones(1,length(brange));
        weights = ones(1,length(brange));
        rinf_mf = zeros(1,length(brange));
        for bi = 1:length(brange)
            X = simRinf{wi,ki,bi};
            [xmean,xsdev,weight] = corrected(X,0.05,0.05);
            weights(bi) = weight/(xsdev+0.1);
            rho_est(bi) = xmean;
            
            rinf_mf(bi)=mfaRinf{wi,ki,bi};
        end
        
        
        % SIMULATION 
        sigmas = 2; % spikeyness of gaussian convolution 
        % (relative importance of the closeby data points, closeness to linear
        % interpolation)
        [ma,bs] = movingaverage(rho_est,weights,brange,betasteps,sigmas);
        atrate_sm(wi,ki) = ma(atindex);
        % MF threshold interpolate
        mf_interp = interp1q(brange',rinf_mf',bs'); 
        atrate_mf(wi,ki) = mf_interp(atindex);
        
    end
end
% load('MFvsSIMin2Datrate.mat') %,'atrate_mf','atrate_sm')

LegendFontsizes = 12;
gamma=1/40;
bc = 0.0025;
wfine=linspace(0.0,0.01,1000);
kapfine=(mu-1)*bc - gamma - wfine;

minval = 0.000;
maxval = 0.6;
figure; 
plot(wrange,zeros(1,length(wrange)));
hold on;
dw=wrange(2)-wrange(1);
dk=kaprange(2)-kaprange(1);
for wi=1:length(wrange)
    w = wrange(wi);
    for ki=1:length(kaprange)
        val = atrate_mf(wi,ki);
        scaled = (val-minval)/(maxval-minval);
        kap = kaprange(ki);
        r = rectangle('Position',[w kap dw dk]');
        r.FaceColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        r.EdgeColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        r.LineWidth = 0.0001;
    end
end

ql=plot(wfine,kapfine);
ql.LineWidth=3;
ql.Color='black';
ql.DisplayName = 'critical curve (MF)';
ax=gca;
ax.YLim=[0,0.01];
ax.XLim=[0,0.01];
ax.XLabel.String='w';
ax.YLabel.String='\kappa';
title('r_{\infty} at \beta=0.0025 (Mean Field)')

colorMap = [linspace(0,1,256)',linspace(0.5,0,256)',linspace(1,0,256)'];
colormap(colorMap);
cb=colorbar;
cb.LineWidth=0.5;
cb.Ticks = [0:0.25:1];
cb.TicksMode ='manual';
cb.TickLabels = minval + [0:0.25:1]*(maxval-minval);

lgd=legend(ql);
lgd.Location='northeast';
lgd.Box='off';
lgd.BoxFace.ColorType='truecoloralpha';
lgd.FontSize=LegendFontsizes;

resolution=300;
folder='figures/';
filename=strcat('MeanField_2d_atrate_date_8Jun');
direction=strcat(folder,filename,'.png');
saveas(gcf,direction)

clear('ql');

figure; 
plot(wrange,zeros(1,length(wrange)));
hold on;
dw=wrange(2)-wrange(1);
dk=kaprange(2)-kaprange(1);
for wi=1:length(wrange)
    w = wrange(wi);
    for ki=1:length(kaprange)
        val = atrate_sm(wi,ki);
        scaled = (val-minval)/(maxval-minval);
        kap = kaprange(ki);
        r = rectangle('Position',[w kap dw dk]');
        r.FaceColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        r.EdgeColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        r.LineWidth = 0.0001;
    end
end
ql=plot(wfine,kapfine);
ql.LineWidth=3;
ql.Color='black';
ql.DisplayName = 'critical curve (MF)';
ax=gca;
ax.YLim=[0,0.01];
ax.XLim=[0,0.01];
ax.XLabel.String='w';
ax.YLabel.String='\kappa';
title('r_{\infty} at \beta=0.0025 (Simulation)')
colorMap = [linspace(0,1,256)',linspace(0.5,0,256)',linspace(1,0,256)'];
colormap(colorMap);
cb=colorbar;
cb.LineWidth=0.5;
cb.Ticks = [0:0.25:1];
cb.TicksMode ='manual';
cb.TickLabels = minval + [0:0.25:1]*(maxval-minval);

lgd=legend(ql);
lgd.Location='northeast';
lgd.Box='off';
lgd.BoxFace.ColorType='truecoloralpha';
lgd.FontSize=LegendFontsizes;

resolution=300;
folder='figures/';
filename=strcat('Simulation_2d_atrate_date_8Jun');

direction=strcat(folder,filename,'.png');
saveas(gcf,direction)