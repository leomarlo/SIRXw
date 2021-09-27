plot_utilities

% load('sweepMay24N500.mat')

%% The big sweep from June
% load('sweepJun06N500delta002.mat')
%% The small sweep from April
load('joinedsweep23AprN500delta002.mat')

atindex = 501;
%% for appeal
atindex = 201;

atrate_mf = zeros(length(wrange),length(kaprange));
atrate_sm = zeros(length(wrange),length(kaprange));
betasteps = 1000;
for wi=1:length(wrange)
    w=wrange(wi)
    for ki=1:length(kaprange)
        kap=kaprange(ki);
        imax_est  = ones(1,length(brange));
        weights = ones(1,length(brange));
        Imax_mf = zeros(1,length(brange));
        for bi = 1:length(brange)
            X = simRinf{wi,ki,bi};
            Imax = simImax{wi,ki,bi};
            [x_mean_X,x_std_X,weight,upperones,condition] = corrected(X,0.05,0.05);
            if condition
                x_mean = mean(Imax(upperones));
                x_std = std(Imax(upperones));
            else
                x_mean = mean(Imax);
                x_std = std(Imax);
            end
            weights(bi) = weight/(x_std+0.1);
            imax_est(bi) = x_mean;
            
            Imax_mf(bi)=mfaImax{wi,ki,bi};
        end
        
        
        % SIMULATION 
        sigmas = 2; % spikeyness of gaussian convolution 
        % (relative importance of the closeby data points, closeness to linear
        % interpolation)
        [imaxs,bs] = movingaverage(imax_est,weights,brange,betasteps,sigmas);
        atrate_sm(wi,ki) = imaxs(atindex);
        % MF threshold interpolate
        mf_interp = interp1q(brange',Imax_mf',bs'); 
        atrate_mf(wi,ki) = mf_interp(atindex);
        
    end
end
% load('MFvsSIMin2Datrate.mat') %,'atrate_mf','atrate_sm')
% load('MFvsSIMin2DatrateForImax.mat') %,'atrate_mf','atrate_sm')
atrate_mf(atrate_mf==0)=nan;
atrate_sm(atrate_sm==0)=nan;

LegendFontsizes = 12;
gamma=1/40;
bc = 0.0025;
wfine=linspace(0.0,0.01,1000);
kapfine=(mu-1)*bc - gamma - wfine;

% for appeal parameters:

%% for the appeal we have
gamma=1/14; 
wfine=linspace(0.0,0.2,1000);
bc = 0.02;
kapfine=(mu-1)*bc - gamma - wfine;


minval = 0.000;
maxval = 0.1;
% minval = 0.000;
minval = round(max(0,min(min(atrate_mf(:)), min(atrate_sm(:))) * ( 1 - 0.05)),2);
maxval = round(min(1,max(max(atrate_mf(:)), max(atrate_sm(:))) * ( 1 + 0.05)),2);


figure; 
plot(wrange,zeros(1,length(wrange)));
hold on;
dw=wrange(2)-wrange(1);
dk=kaprange(2)-kaprange(1);
for wi=1:length(wrange)
    w = wrange(wi);
    for ki=1:length(kaprange)
        val = atrate_mf(wi,ki);
        kap = kaprange(ki);
        r = rectangle('Position',[w kap dw dk]');
        if isnan(val)
            r.FaceColor = 'white';
            r.EdgeColor = 'white';
        else
            scaled = (val-minval)/(maxval-minval);
            r.FaceColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
            r.EdgeColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        end
        r.LineWidth = 0.0001;
    end
end

ql=plot(wfine,kapfine);
ql.LineWidth=3;
ql.Color='black';
ql.DisplayName = 'critical curve (MF)';
ax=gca;
ax.YLim=[0,max(kaprange)];
ax.XLim=[0,max(wrange)];
ax.XLabel.String='w';
ax.YLabel.String='\kappa';
ax.XTick = ticks2d.*20;
ax.YTick = ticks2d.*20;
ax.TickDir = 'out';
ax.FontSize = fontsize;
ax.LineWidth = FrameThickness;
ax.TickLength=TickLengths;
ax.FontWeight = FontWeights;
title({'Mean Field'})

set(ax,'box','off');
grid on
daspect([1 1 1])


% title('r_{\infty} at \beta=0.0025 (Mean Field)')
colorMap = [linspace(0,1,256)',linspace(0.5,0,256)',linspace(1,0,256)'];
colormap(colorMap);
cb=colorbar;
cb.LineWidth=0.5;
cb.Ticks = [0:0.25:1];
cb.TicksMode ='manual';
cb.TickLabels = minval + [0:0.25:1]*(maxval-minval);
title(cb, 'r_\infty');
cb.FontSize = fontsize; 
cb.LineWidth = FrameThickness;
cb.FontWeight = FontWeights;
pos=get(ax,'position');  % retrieve the current values
pos(4)=0.97*pos(4);        % try reducing width 10%
pos(3)=0.95*pos(3);  
set(ax,'position',pos);

lgd=legend(ql);
lgd.Location='northeast';
lgd.Box='off';
lgd.BoxFace.ColorType='truecoloralpha';
lgd.FontSize=LegendFontsize;

resolution=300;
folder='figures/';
% filename=strcat('MeanField_2d_atrate_delta2_Imax_date_18Jun');
filename=strcat('atrate_2d_Imaximal_MF_',date);
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
        kap = kaprange(ki);
        r = rectangle('Position',[w kap dw dk]');
        if isnan(val)
            r.FaceColor = 'white';
            r.EdgeColor = 'white';
        else
            scaled = (val-minval)/(maxval-minval);
            r.FaceColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
            r.EdgeColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        end
        r.LineWidth = 0.0001;
    end
end
ql=plot(wfine,kapfine);
ql.LineWidth=3;
ql.Color='black';
ql.DisplayName = 'critical curve (MF)';
ax=gca;
ax.YLim=[0,max(kaprange)];
ax.XLim=[0,max(wrange)];
ax.XLabel.String='w';
ax.YLabel.String='\kappa';
ax.XTick = ticks2d.*20;
ax.YTick = ticks2d.*20;
ax.TickDir = 'out';
ax.FontSize = fontsize;
ax.LineWidth = FrameThickness;
ax.TickLength=TickLengths;
ax.FontWeight = FontWeights;
title({'Simulation'})

set(ax,'box','off');
grid on
daspect([1 1 1])

% title('r_{\infty} at \beta=0.0025 (Simulation)')

pos=get(ax,'position');  % retrieve the current values
pos(4)=0.97*pos(4);        % try reducing width 10%
pos(3)=0.95*pos(3);  
set(ax,'position',pos);

lgd=legend(ql);
lgd.Location='northeast';
lgd.Box='off';
lgd.BoxFace.ColorType='truecoloralpha';
lgd.FontSize=LegendFontsize;

resolution=300;
folder='figures/';
% filename=strcat('Simulation_2d_atrate_delta2_Imax_date_18Jun');

filename=strcat('atrate_2d_Imaximal_SIM_',date);

direction=strcat(folder,filename,'.png');
saveas(gcf,direction)