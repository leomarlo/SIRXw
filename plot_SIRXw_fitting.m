% load('sweepMay24N500.mat')

critical=zeros(length(wrange),length(kaprange));

mfImaxatbeta=zeros(length(wrange),length(kaprange));
mfRinfatbeta=zeros(length(wrange),length(kaprange));

nmImaxatbeta=zeros(length(wrange),length(kaprange));
nmRinfatbeta=zeros(length(wrange),length(kaprange));

betasteps = 1000;
beta_fine = brange(1):(brange(end)-brange(1))/(betasteps-1):brange(end);

% weighted_moving_average = zeros(length(wrange),length(kaprange), length(beta_fine));
 
% for wi=1:length(wrange)
%     w=wrange(wi)
%     for ki=1:length(kaprange)
%         kap=kaprange(ki);            
%             for bi = 1:length(brange)
%                 beta = brange(bi);
%                 
%             end
%     end
% end

sample = [[2,2];[4,4];[13,25];[12,51];[51,12];[51,51]];

% sample = [[2,2];[3,3]]

for s=sample'
    wi = s(1); w = wrange(wi)
    ki = s(2); kap = kaprange(ki) 


    mean_rinf = zeros(1,length(brange));
    mean_rinf_mf = zeros(1,length(brange));
    sdev_rinf = zeros(1,length(brange));
    for bi = 1:length(brange)
        mean_rinf(bi)=mean(simRinf{wi,ki,bi});
        sdev_rinf(bi)=std(simRinf{wi,ki,bi}); 
        mean_rinf_mf(bi)=mfaRinf{wi,ki,bi};
    end
    % pl(1) = errorbar(brange,mean_rinf,sdev_rinf);
    % pl(1).Marker = 'o';
    % pl(1).MarkerSize = 2; pl(1).LineWidth = 2;
    % hold on;
    % pl(2) = plot(brange,mean_rinf_mf);
    % pl(2).LineWidth = 2;
    % 
    % clear('pl')

    % weighted_moving_average = ones(1,length(beta_fine));

    rinfs = zeros(length(brange),length(simRinf{wi,ki,1}));
    plot_bool = zeros(length(brange),length(simRinf{wi,ki,1}));
    figure;
    pl(1) = plot(brange,mean_rinf_mf);
    pl(1).LineWidth = 3;
    pl(1).DisplayName = 'Mean Field';
    % pl(1).LineStyle = '--';
    pl(1).Color = 'black';
    hold on;
    pl(6) = plot(brange, 0.05 * ones(1,length(brange)));
    pl(6).LineStyle = '--';
    pl(6).LineWidth = 2;
    pl(6).Color = [0.9290 0.6940 0.1250];
    pl(6).DisplayName = '\tau_{Kurtosis}';
    dbeta = brange(2)-brange(1);
    rho_est  = ones(1,length(brange));
    weights = ones(1,length(brange));
    for bi = 1:length(brange)
        bla = simRinf{wi,ki,bi};
        if std(bla) <0.000001
            idx = ones(length(bla),1);
            C=[bla(1),bla(1)];
            Ku =0;
%             [idx,C] = kmeans(bla',2);
        else
            [idx,C] = kmeans(bla',2);
            Ku = kurtosis(bla)*(std(bla)^2);
        end
        
    %     JB = JarqueBera(bla);
    %     plot(brange(bi),JB/10,'*','MarkerFaceColor','blue','MarkerEdgeColor','blue','MarkerSize',4);   
        
        pl(2) = plot(brange(bi),Ku*3,'sq','MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerSize',7,'LineWidth',4);
        pl(2).DisplayName = 'Unnormalized Kurtosis';
        if C(2) < C(1)
            idx = 3-idx;
            pl(3) = plot(brange(bi),C(1),'*','MarkerEdgeColor','blue','MarkerSize',6,'LineWidth',4);
            pl(4) = plot(brange(bi),C(2),'*','MarkerEdgeColor','green','MarkerSize',6,'LineWidth',4);
        else
            pl(4) = plot(brange(bi),C(1),'*','MarkerEdgeColor','green','MarkerSize',6,'LineWidth',4);
            pl(3) = plot(brange(bi),C(2),'*','MarkerEdgeColor','blue','MarkerSize',6,'LineWidth',4);
        end
        pl(3).DisplayName = '2-means upper centroid';
        pl(4).DisplayName = '2-means lower centroid';
        [xmean,xsdev,weight] = corrected(bla,0.05,0.05);
        weights(bi) = weight/(xsdev+0.1);
        rho_est(bi) = xmean;
        errbarwidth = 2;
        errbarcol = 'red';
        pl(5)  = plot (brange(bi),xmean);
        pl(5).Marker ='o';
        pl(5).MarkerSize = 4;
        pl(5).MarkerEdgeColor = errbarcol;
        pl(5).LineWidth = 3;
        pl(5).LineStyle = 'none';
        pl(5).DisplayName = 'surviving runs';
        bars0  = plot ([brange(bi),brange(bi)],[xmean+xsdev,xmean-xsdev]);
        bars0.LineStyle = '-';
        bars0.LineWidth = errbarwidth;
        bars0.Color =errbarcol;
        bars1  = plot ([brange(bi)-dbeta/8,brange(bi)+dbeta/8],[xmean+xsdev,xmean+xsdev]);
        bars1.LineStyle = '-';
        bars1.LineWidth = errbarwidth;
        bars1.Color =errbarcol;
        bars2  = plot ([brange(bi)-dbeta/8,brange(bi)+dbeta/8],[xmean-xsdev,xmean-xsdev]);
        bars2.LineStyle = '-';
        bars2.LineWidth = errbarwidth;
        bars2.Color =errbarcol;
        for sa = 1:length(simRinf{wi,ki,1})
            rinfs(bi,sa) =  bla(sa);

            plot_bool(bi,sa) = idx(sa);
            if idx(sa)==1
                plot(brange(bi),rinfs(bi,sa),'.','MarkerFaceColor','green','MarkerEdgeColor','green')
            else
                plot(brange(bi),rinfs(bi,sa),'.','MarkerFaceColor','blue','MarkerEdgeColor','blue')
            end
        end
    end
    sigmas = 2; % spikeyness of gaussian convolution 
    % (relative importance of the closeby data points, closeness to linear
    % interpolation)
    [ma,bs] = movingaverage(rho_est,weights,brange,betasteps,sigmas);
    pl(7) = plot(bs,ma);
    pl(7).DisplayName = 'Weighted Moving Average';
    pl(7).Color = [0.5,0.5,0.5];
    pl(7).LineWidth = 4;

    ax = gca;
    ax.YLim = [0,1];
    ax.XLim = [-0.0001,0.0051];
    ax.XLabel.String = '\beta';
    ax.YLabel.String = 'r_\infty';
    titletext = strcat('Simulation vs MF (kappa= ',num2str(kap),' and w= ',num2str(w),')');
    title(titletext);


    lgd = legend(pl);
    lgd.Location='northwest';
    
    resolution=300;
    folder='figures/';
    filename=strcat('MF_vs_sim_kap_',num2str(kap),'_w_',num2str(w),'_date_8Jun');

    direction=strcat(folder,filename,'.png');
    saveas(gcf,direction)

end


% % figure;
% mean_rinf = zeros(1,length(brange));
% mean_rinf_mf = zeros(1,length(brange));
% sdev_rinf = zeros(1,length(brange));
% for bi = 1:length(brange)
%     mean_rinf(bi)=mean(simRinf{wi,ki,bi});
%     sdev_rinf(bi)=std(simRinf{wi,ki,bi}); 
%     mean_rinf_mf(bi)=mfaRinf{wi,ki,bi};
% end
% % pl(1) = errorbar(brange,mean_rinf,sdev_rinf);
% % pl(1).Marker = 'o';
% % pl(1).MarkerSize = 2; pl(1).LineWidth = 2;
% % hold on;
% % pl(2) = plot(brange,mean_rinf_mf);
% % pl(2).LineWidth = 2;
% % 
% % clear('pl')
% 
% % weighted_moving_average = ones(1,length(beta_fine));
% 
% rinfs = zeros(length(brange),length(simRinf{wi,ki,1}));
% plot_bool = zeros(length(brange),length(simRinf{wi,ki,1}));
% figure;
% pl(1) = plot(brange,mean_rinf_mf);
% pl(1).LineWidth = 3;
% pl(1).DisplayName = 'Mean Field';
% % pl(1).LineStyle = '--';
% pl(1).Color = 'black';
% hold on;
% pl(6) = plot(brange, 0.05 * ones(1,length(brange)));
% pl(6).LineStyle = '--';
% pl(6).LineWidth = 2;
% pl(6).Color = [0.9290 0.6940 0.1250];
% pl(6).DisplayName = '\tau_{Kurtosis}';
% dbeta = brange(2)-brange(1);
% rho_est  = ones(1,length(brange));
% weights = ones(1,length(brange));
% for bi = 1:length(brange)
%     bla = simRinf{wi,ki,bi};
%     [idx,C] = kmeans(bla',2);
% %     JB = JarqueBera(bla);
% %     plot(brange(bi),JB/10,'*','MarkerFaceColor','blue','MarkerEdgeColor','blue','MarkerSize',4);   
%     Ku = kurtosis(bla)*(std(bla)^2);
%     pl(2) = plot(brange(bi),Ku*3,'sq','MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerSize',7,'LineWidth',4);
%     pl(2).DisplayName = 'Unnormalized Kurtosis';
%     if C(2) < C(1)
%         idx = 3-idx;
%         pl(3) = plot(brange(bi),C(1),'*','MarkerEdgeColor','blue','MarkerSize',6,'LineWidth',4);
%         pl(4) = plot(brange(bi),C(2),'*','MarkerEdgeColor','green','MarkerSize',6,'LineWidth',4);
%     else
%         pl(4) = plot(brange(bi),C(1),'*','MarkerEdgeColor','green','MarkerSize',6,'LineWidth',4);
%         pl(3) = plot(brange(bi),C(2),'*','MarkerEdgeColor','blue','MarkerSize',6,'LineWidth',4);
%     end
%     pl(3).DisplayName = '2-means upper centroid';
%     pl(4).DisplayName = '2-means lower centroid';
%     [xmean,xsdev,weight] = corrected(bla,0.05,0.05);
%     weights(bi) = weight/xsdev;
%     rho_est(bi) = xmean;
%     errbarwidth = 2;
%     errbarcol = 'red';
%     pl(5)  = plot (brange(bi),xmean);
%     pl(5).Marker ='o';
%     pl(5).MarkerSize = 4;
%     pl(5).MarkerEdgeColor = errbarcol;
%     pl(5).LineWidth = 3;
%     pl(5).LineStyle = 'none';
%     pl(5).DisplayName = 'surviving runs';
%     bars0  = plot ([brange(bi),brange(bi)],[xmean+xsdev,xmean-xsdev]);
%     bars0.LineStyle = '-';
%     bars0.LineWidth = errbarwidth;
%     bars0.Color =errbarcol;
%     bars1  = plot ([brange(bi)-dbeta/8,brange(bi)+dbeta/8],[xmean+xsdev,xmean+xsdev]);
%     bars1.LineStyle = '-';
%     bars1.LineWidth = errbarwidth;
%     bars1.Color =errbarcol;
%     bars2  = plot ([brange(bi)-dbeta/8,brange(bi)+dbeta/8],[xmean-xsdev,xmean-xsdev]);
%     bars2.LineStyle = '-';
%     bars2.LineWidth = errbarwidth;
%     bars2.Color =errbarcol;
%     for sa = 1:length(simRinf{wi,ki,1})
%         rinfs(bi,sa) =  bla(sa);
%         
%         plot_bool(bi,sa) = idx(sa);
%         if idx(sa)==1
%             plot(brange(bi),rinfs(bi,sa),'.','MarkerFaceColor','green','MarkerEdgeColor','green')
%         else
%             plot(brange(bi),rinfs(bi,sa),'.','MarkerFaceColor','blue','MarkerEdgeColor','blue')
%         end
%     end
% end
% sigmas = 2; % spikeyness of gaussian convolution 
% % (relative importance of the closeby data points, closeness to linear
% % interpolation)
% [ma,bs] = movingaverage(rho_est,weights,brange,betasteps,sigmas);
% pl(7) = plot(bs,ma);
% pl(7).DisplayName = 'Weighted Moving Average';
% pl(7).Color = [0.5,0.5,0.5];
% pl(7).LineWidth = 4;
% 
% ax = gca;
% ax.YLim = [0,1];
% ax.XLim = [-0.0001,0.0051];
% ax.XLabel.String = '\beta';
% ax.YLabel.String = 'r_\infty';
% titletext = strcat('Simulation vs MF (kappa= ',num2str(kap),' and w= ',num2str(w),')');
% title(titletext);
% 
% 
% lgd = legend(pl);
% lgd.Location='northwest';
% 
% figure;
% for sa = 1:length(simRinf{wi,ki,1})
%     q(sa) = plot(brange, rinfs(:,sa));
%     q(sa).LineStyle = 'none';
%     q(sa).Marker = 'o';
%     q(sa).MarkerSize = 2;
%     q(sa).LineWidth = 3;
%     hold on;
% end
% q(sa+1) = plot(brange,mean_rinf_mf);
% q(sa+1).LineWidth = 2;
% 
%  lgd=legend(q);
%  lgd.Location='northwest';
%  
% resolution=300;
% folder='figures/';
% filename=strcat('Jun08_MF_vs_simulation');
% 
% direction=strcat(folder,filename,'.png');
% saveas(gcf,direction)
