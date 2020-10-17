tmax = 1000;
gamma = 1/40;
delta = 1/10;
mu = 15;
wrange = 0:0.002:0.035;
krange = 0:0.0001:0.035;
brange = [0.0020,0.0025,0.003];
effconn = [-0.8,0,5];
xlim = [min(wrange),max(wrange)];
ylim = [min(krange),max(krange)];
% 
% rho0 = 0.001;
% rhosi0 = rho0 * mu *(1+ eff) ;
% rhoss0 = (mu/2)- rhosi0;
% ini=[1-rho0,rho0,0,rhosi0,rhoss0];

betacrit = @(w,k,gamma,mu) (w+k+gamma)/(mu+1);
kcrits = @(beta,gamma,w,mu) (beta*(mu-1) - gamma - w);

non_poss = ones(length(wrange),length(effconn),length(brange));
crit_mfs = ones(length(wrange),length(effconn),length(brange));

for ei = 1:length(effconn)
    eff = effconn(ei);
    rho0 = 0.001;
    rhosi0 = rho0 * mu *(1+ eff) * (1-rho0);
    rhoss0 = (mu/2)- rhosi0;
    ini=[1-rho0,rho0,0,rhosi0,rhoss0];
    for bi = 1:length(brange)
        beta = brange(bi);
        disp(beta);
        non_pos = ones(1,length(wrange));
        crit_mf = ones(1,length(wrange));
        parfor wi = 1:length(wrange)
            
            w = wrange(wi);
            non_positive = NaN;
            for ki = 1:length(krange)
                k = krange(ki);
                ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,k,delta);
                [ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);

                fsi = -(beta+gamma+w+k)*xs(:,4) + beta*xs(:,4).*(2*xs(:,5)-xs(:,4))./xs(:,1);
                fss = -2*beta*xs(:,4).*xs(:,5)./xs(:,1) + w*(xs(:,1)./(xs(:,1)+xs(:,3))).*xs(:,4);
                diff = fsi - 2*fss;
                if any(diff<0.000000000001)
                    if ki == 1
                        non_positive = NaN;
                    else
                        non_positive = k;
                    end
                    break
                end
            end
            non_pos(wi) = non_positive;
            crit_mf(wi) = kcrits(beta,gamma,w,mu);
        end
        non_poss(:,ei,bi)=non_pos;
        crit_mfs(:,ei,bi)=crit_mf;
    end
end

% 
% 
% cols = rand(length(brange),3);
% styles = {':','--','-'};
% cols = eye(3);
% pl = 1;
% 
% figure;
% for ei = 1:length(effconn)
%     clear('pl')
%     for bi = 1:length(brange)
%         
%         pl(2*(bi-1)+1) = plot(wrange,non_poss(:,ei,bi));
%         pl(2*(bi-1)+1).Color = cols(bi,:);
%         pl(2*(bi-1)+1).LineStyle = styles{ei};
%         pl(2*(bi-1)+1).DisplayName = strcat('non-negative (\beta = ',num2str(brange(bi)),')');
%         pl(2*(bi-1)+1).LineWidth = 2.5;
%         hold on;
%         pl(2*bi) = plot(wrange,crit_mfs(:,ei,bi));
%         pl(2*bi).Color = cols(bi,:);
%         pl(2*bi).DisplayName = strcat('MF transition (\beta = ',num2str(brange(bi)),')');
%         pl(2*bi).LineStyle = '-';
%         pl(2*bi).LineWidth = 1.5;
%     end
%     ax = gca; 
%     ax.XLim = xlim;
%     ax.YLim = ylim;
% end
% 
% lgd = legend(pl);
% lgd.Location = 'northeast';
% resolution=300;
% 
% folder='figures/';
% filename=strcat('positivity_23Aug');
% 
% direction=strcat(folder,filename,'.png');
% saveas(gcf,direction)
% 


cols = rand(length(brange),3);
styles = {':','--','-.'};
% cols = eye(3);
cols = [[1,0,0];[0,0,0];[0.5,0.5,0.5];[0,0,1]];
pl = 1;

figure;
clear('pl')
for ei = 1:length(effconn)
    
    for bi = 1:length(brange)
        
        pl(ei) = plot(wrange,non_poss(:,ei,bi));
        pl(ei).Color = cols(bi,:);
        pl(ei).LineStyle = styles{ei};
        pl(ei).DisplayName = strcat('positivity @ \rho_{SI}(0) = ',num2str((1+effconn(ei))*mu),'\rho_{S}(0)');
        pl(ei).LineWidth = 2.5;
        hold on;
        pl(length(effconn)+(length(brange)-bi+1)) = plot(wrange,crit_mfs(:,ei,bi));
        pl(length(effconn)+(length(brange)-bi+1)).Color = cols(bi,:);
        pl(length(effconn)+(length(brange)-bi+1)).DisplayName = strcat('MF transition @ \beta=',num2str(brange(bi)),'');
        pl(length(effconn)+(length(brange)-bi+1)).LineStyle = '-';
        pl(length(effconn)+(length(brange)-bi+1)).LineWidth = 1.5;
    end
    
    ppll = plot(wrange,non_poss(:,ei,bi));
    ppll.Color = cols(bi+1,:);
    ppll.LineStyle = styles{ei};
    ppll.DisplayName = strcat('positivity\rho_{SI}(0) = ',num2str((1+effconn(ei))*mu),'\rho_{S}(0)');
    ppll.LineWidth = 2.5;
%     qqll = plot(wrange,crit_mfs(:,ei,bi));
%     qqll.Color = cols(bi+1,:);
%     qqll.DisplayName = strcat('MF transition (\beta=',num2str(brange(bi)),')');
%     qqll.LineStyle = '-';
%     qqll.LineWidth = 1.5;
    pl(length(effconn)+(length(brange)-bi+1)) = plot(wrange,crit_mfs(:,ei,bi));
    pl(length(effconn)+(length(brange)-bi+1)).Color = cols(bi+1,:);
    pl(length(effconn)+(length(brange)-bi+1)).DisplayName = strcat('MF transition @ \beta=',num2str(brange(bi)),'');
    pl(length(effconn)+(length(brange)-bi+1)).LineStyle = '-';
    pl(length(effconn)+(length(brange)-bi+1)).LineWidth = 1.5;
    ax = gca; 
    ax.XLim = xlim;
    ax.YLim = ylim;
end


lgd = legend(pl);
% lgd = legend(pl(1:(length(effconn)+1)));
lgd.Location = 'northeast';
resolution=300;

folder='figures/';
filename=strcat('positivity_10Sep');

direction=strcat(folder,filename,'.png');
saveas(gcf,direction)

