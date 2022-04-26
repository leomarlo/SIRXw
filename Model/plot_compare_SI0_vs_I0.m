% COMPARE R INFTY AS FUNCTION OF I0 AND SI0

%%% load the data from test run %%%%
load('i0svssi0s.mat')

%%% load parameters 
% N=1000; mu=15; L = ceil(mu*N/2);
tmax = 1500;
beta=0.004; gamma=1/40; w=1/100; 
kap=12/10000; delta=2/100;
betarange = [0.002:0.0001:0.006];
% rho0=1/100;
SA = 100; 
Nt = 1000;

% initialize plots
pl = [];

%% COMPARISONS
%% compare low and high si0 at the same i0
casenumber = 1; % low i0 low si0
casenumber = 3; % low i0 high si0


means1 = zeros(2, length(betarange),SA);
means2 = zeros(2, length(betarange),SA);
mf = zeros(2, length(betarange));
theory1 = zeros(2, length(betarange));
theory2 = zeros(2, length(betarange));

cc = 1;
for casenumber = [1,3]
    disp('start first case')
    disp(res{casenumber}.what)
    %% load A and x
    A = res{casenumber}.A;
    x = res{casenumber}.x;
    B = res{casenumber}.B;
    y = res{casenumber}.y;
    N = size(A,1);

    %% compute initial values for ODE integration
    si0=((x==0)*A)*(x==1)'/N;
    i0 =sum(x)/N;
    rhos=1-i0; rhor=0;
    rhoss= (mu/2)-si0;
    ini=[rhos,i0,rhor,si0,rhoss];
    disp(strcat('the computed si0 is: ',num2str(si0),' and it should be ',num2str(res{casenumber}.si0)))
    disp(strcat('the computed si0 is: ',num2str(i0),' and it should be ',num2str(res{casenumber}.i0)))

    %% ODE integration
    ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,kap,delta);
    [ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);
    Tt= min(Nt,length(ts));

    %% compute result from A and x
    bi = 1;
    for beta=betarange
        disp('beta')
        disp(beta)
        rinf1 = zeros(1, SA);
        rinf2 = zeros(1, SA);
        ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,kap,delta);
        [ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);
        Tt= min(Nt,length(ts));
        for sa=1:SA
%             disp(sa)
            result1 = SIRXi_w_fromAandx(A,x,beta,gamma,w,kap,delta,tmax,Nt); 
            result2 = SIRXi_w_fromAandx(B,y,beta,gamma,w,kap,delta,tmax,Nt); 
            rinf1(sa) = (result1.NRs(end) + result1.NXis(end))/result1.N;
            rinf2(sa) = (result2.NRs(end) + result2.NXis(end))/result2.N;
        end

        means1(cc,bi,:) = rinf1;
        means2(cc,bi,:) = rinf2;
        mf(cc,bi) = xs(end,3);
        theory1(cc,bi) = i0 + (beta/(beta + gamma + w + kap))*si0;
        theory2(cc,bi) = (gamma/(gamma + kap))*i0 + ((gamma*beta^2)/((gamma+kap)*(beta + gamma + w + kap)))*si0;
        bi = bi +1;
    end
    cc = cc +1;
end
% beta
%% plot 

casenumbers = [1,3];
for cc=1:2
    casenumber = casenumbers(cc);
    figure; 
    bi = 1;
    for beta=betarange
        at_beta_1 = reshape(means1(cc,bi,:),[1,SA]);
        po1 = plot( beta*ones(1,length(at_beta_1)),at_beta_1); 
        po1.Marker = '.';
        po1.LineStyle = 'none';
        po1.Color = [0,0,0];
        po1.DisplayName = 'Simulation 1';
        hold on;
        at_beta_2 = reshape(means2(cc,bi,:),[1,SA]);
        po2 = plot( beta*ones(1,length(at_beta_2)),at_beta_2); 
        po2.Marker = '.';
        po2.LineStyle = 'none';
        po2.Color = [0.5,0.5,0.5];
        po2.DisplayName = 'Simulation 2';
        bi = bi+1;
    end
    po3 = plot(betarange, mf(cc,:));
    po3.Color = [1,0,0.5];
    po3.LineStyle = '-';
    po3.LineWidth = 2;
    po3.DisplayName = 'Mean Field';
    po4 = plot(betarange, theory1(cc,:));
    po4.Color = [0.5,0,1];
    po4.LineStyle = ':';
    po4.LineWidth = 2;
    po4.DisplayName = 'Theory 1';
    po5 = plot(betarange, theory2(cc,:));
    po5.Color = [0.5,1,0];
    po5.LineStyle = ':';
    po5.LineWidth = 2;
    po5.DisplayName = 'Theory 2';

    ax = gca;
    ax.YLabel.String = '\rho_R(\infty)';
    ax.XLabel.String = '\beta';
    title(strcat('Theory vs Simulation (\rho_I(0)=',num2str(res{casenumber}.i0),', \rho_{SI}(0)=',num2str(res{casenumber}.si0),')'));

    lgd = legend([po1,po2,po3,po4,po5]);
    lgd.Location = 'northwest';

    resolution=600;
    folder='figures/';
    filename=strcat('Jul05_MF_vs_simulation_vs_theory_',res{casenumber}.what);

    direction=strcat(folder,filename,'.png');
    saveas(gcf,direction)
end

%% both

figure; 
bi = 1;
for beta=betarange
    cc=1;
    at_beta_1 = reshape(means1(1,bi,:),[1,SA]);
    plo1 = plot( beta*ones(1,length(at_beta_1)),at_beta_1); 
    plo1.Marker = 'o';
    plo1.MarkerSize = 1;
    plo1.LineStyle = 'none';
    plo1.Color = [0.85,0,0];
    plo1.DisplayName = strcat('Simulation (\rho_{SI}(0)=',num2str(res{casenumbers(cc)}.si0),')');
    hold on;
    at_beta_2 = reshape(means2(1,bi,:),[1,SA]);
    plo2 = plot( beta*ones(1,length(at_beta_2)),at_beta_2); 
    plo2.Marker = 'o';
    plo2.MarkerSize = 1;
    plo2.LineStyle = 'none';
    plo2.Color = [0.85,0.0,0.0];
    plo2.DisplayName = 'Simulation';

    cc=2;
    at_beta_1 = reshape(means1(2,bi,:),[1,SA]);
    plo3 = plot( beta*ones(1,length(at_beta_1)),at_beta_1); 
    plo3.Marker = '.';
    plo3.MarkerSize = 3;
    plo3.LineStyle = 'none';
    plo3.Color = [0,0,0.85];
    plo3.DisplayName = strcat('Simulation (\rho_{SI}(0)=',num2str(res{casenumbers(cc)}.si0),')');
    at_beta_2 = reshape(means2(2,bi,:),[1,SA]);
    plo4 = plot( beta*ones(1,length(at_beta_2)),at_beta_2); 
    plo4.Marker = '.';
    plo4.MarkerSize = 3;
    plo4.LineStyle = 'none';
    plo4.Color = [0.0,0.0,0.85];
    plo4.DisplayName = 'Simulation';
    bi = bi+1;
end
cc=1;
plo5 = plot(betarange, mf(cc,:));
plo5.Color = [1,0,0.2];
plo5.LineStyle = '-';
plo5.LineWidth = 2;
plo5.DisplayName = strcat('Mean Field (\rho_{SI}(0)=',num2str(res{casenumbers(cc)}.si0),')');

cc =2;
plo6 = plot(betarange, mf(cc,:));
plo6.Color = [0.2,0,1];
plo6.LineStyle = '-';
plo6.LineWidth = 2;
plo6.DisplayName = strcat('Mean Field (\rho_{SI}(0)=',num2str(res{casenumbers(cc)}.si0),')');

ax = gca;
ax.YLabel.String = '\rho_R(\infty)';
ax.XLabel.String = '\beta';
title(strcat('Theory vs Simulation (low and high initial SI-link-density)'));

lgd = legend([plo1,plo3,plo5,plo6]);
lgd.Location = 'east';

resolution=600;
folder='figures/';
filename=strcat('Jul05_MF_vs_simulation_two_si0_comparison');

direction=strcat(folder,filename,'.png');
saveas(gcf,direction)

%%
% no theory plots 

casenumbers = [1,3];
for cc=1:2
    casenumber = casenumbers(cc);
    figure; 
    bi = 1;
    for beta=betarange
        at_beta_1 = reshape(means1(cc,bi,:),[1,SA]);
        po1 = plot( beta*ones(1,length(at_beta_1)),at_beta_1); 
        po1.Marker = '.';
        po1.LineStyle = 'none';
        po1.Color = [0,0,0];
        po1.DisplayName = 'Simulation 1';
        hold on;
        at_beta_2 = reshape(means2(cc,bi,:),[1,SA]);
        po2 = plot( beta*ones(1,length(at_beta_2)),at_beta_2); 
        po2.Marker = '.';
        po2.LineStyle = 'none';
        po2.Color = [0.5,0.5,0.5];
        po2.DisplayName = 'Simulation 2';
        bi = bi+1;
    end
    po3 = plot(betarange, mf(cc,:));
    po3.Color = [1,0,0.5];
    po3.LineStyle = '-';
    po3.LineWidth = 2;
    po3.DisplayName = 'Mean Field';

    ax = gca;
    ax.YLabel.String = '\rho_R(\infty)';
    ax.XLabel.String = '\beta';
    title(strcat('Theory vs Simulation (\rho_I(0)=',num2str(res{casenumber}.i0),', \rho_{SI}(0)=',num2str(res{casenumber}.si0),')'));

    lgd = legend([po1,po2,po3]);
    lgd.Location = 'northwest';

    resolution=600;
    folder='figures/';
    filename=strcat('Jul05_MF_vs_simulation_vs_',res{casenumber}.what);

    direction=strcat(folder,filename,'.png');
    saveas(gcf,direction)
end

%% compare low and high i0 at the same si0
