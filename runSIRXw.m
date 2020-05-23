tmax=10000;
Nt=1000;

N=1000;
mu=6;
beta=0.001;
gamma=0.002;
w=0.0005;
rho0=0.02;
kap0=0.0001;
kap=0.0004;

% result=SIRXw(N,mu,beta,gamma,w,kap0,kap,rho0,tmax,Nt);
% 
% pl(1)=plot(result.times, result.NIs/result.N);pl(1).Color='red';
% hold on;
% pl(2)=plot(result.times, result.NRs/result.N);pl(2).Color='blue';
% pl(3)=plot(result.times, result.NSs/result.N);pl(3).Color='green';
% pl(4)=plot(result.times, result.NXs/result.N);pl(4).Color='black';
% pl(5)=plot(result.times, (result.NSs+result.NIs+result.NRs+result.NXs)/result.N);
% ax=gca;
% ax.YLim=[0,1];
% lgd=legend(pl);
% lgd.Location='northeast';




% lgd.Orientation='vertical';
% lgd.Box='off';
% lgd.FontSize=LegendFontsizes*(0.5);
% lgd.Color=blackcolor;
% lgd.LineWidth=LineWidths;

%% RUN A PARAMETER SWEEP

wrange=0:0.00005:0.003;
K=10;
    maxI=zeros(length(wrange),K);
    limR=zeros(length(wrange),K);
    for wind=1:length(wrange)
        w=wrange(wind);
        disp(w)
        for sa=1:K
            disp(sa)
            result=SIRXw(N,mu,beta,gamma,w,kap0,kap,rho0,tmax,Nt);
            maxI(wind,sa)=max(result.NIs/N);
            limR(wind,sa)=max(result.NRs/N);

        end
    end

    sweep.maxI=maxI;
    sweep.Rinf=limR;

    folder='figures/';
    filename=strcat('Rinf_w_sweep_16May_Model_','');
    resol=300;
    resolution=strcat('-r',num2str(resol));

    xs=sweep.Rinf;
    fig=figure;
    er=errorbar(wrange,mean(xs,2),std(xs'));
    er.Marker='o';
    er.MarkerSize=6;
    ax=gca;
    ax.YLim=[0,1];
    ax.XLim=[0,max(wrange)];
    ax.XLabel.String='w';
    ax.YLabel.String='r_{\infty}';

    direction=strcat(folder,filename,'.png');
    print(direction,'-dpng',resolution);
%     savefig(fig,strcat(direction(1:(end-4)),'.fig'));

    folder='figures/';
    filename=strcat('maxI_w_sweep_16May_Model_','');
    resol=300;
    resolution=strcat('-r',num2str(resol));

    xs=sweep.maxI;
    fig=figure;
    er=errorbar(wrange,mean(xs,2),std(xs'));
    er.Marker='o';
    er.MarkerSize=6;
    ax=gca;
    ax.YLim=[0,0.5];
    ax.XLim=[0,max(wrange)];
    ax.XLabel.String='w';
    ax.YLabel.String='max(\rho_I)';

    direction=strcat(folder,filename,'.png');
    print(direction,'-dpng',resolution);
%     savefig(fig,strcat(direction(1:(end-4)),'.fig'));

    save(strcat('SIRXw_16May.mat'),'sweep')

