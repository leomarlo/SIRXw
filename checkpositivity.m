
tmax = 1000;
gamma = 1/40;
delta = 1/10;
mu = 15;
wrange = 0:0.005:0.05;
krange = 0:0.005:0.05;
brange = 0.001:0.03;
rho0 = 0.001;
rhosi0 = rho0 * mu * (1-rho0);
rhoss0 = (mu/2)- rhosi0;
ini=[1-rho0,rho0,0,rhosi0,rhoss0];

betacrit = @(w,k,gamma,mu) (w+k+gamma)/(mu+1);
betacrit(0.05,0.05,gamma,mu)

valid = ones(length(wrange),length(krange),length(brange));
for wi = 1:length(wrange)
    w = wrange(wi);
    disp('rewiring is')
    disp(w)
    for ki = 1:length(wrange)
        k = krange(ki);
        valid_bs = ones(1,length(brange));
        parfor bi = 1:length(brange)
            beta = brange(bi);
            ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,k,delta);
            [ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);

            fsi = -(beta+gamma+w+k)*xs(:,4) + beta*xs(:,4).*(2*xs(:,5)-xs(:,4))./xs(:,1);
            fss = -2*beta*xs(:,4).*xs(:,5)./xs(:,1) + w*(xs(:,1)./(xs(:,1)+xs(:,3))).*xs(:,4);
            
            diff = fsi - 2*fss;
            sdiff = sum(diff);
            if sdiff < 0
                valid_bs(bi) = - sum(diff(diff>0))/sum(abs(diff));
            else
                valid_bs(bi) = - sum(diff(diff<0))/sum(abs(diff));
            end
        end
        valid(wi,ki,:) = valid_bs;
        
        if mod(ki,5)==0
            
            disp('quarantine is')
            disp(k)
       
        end
    end
end
% 
% 
ki = 10; k = krange(ki);
wi = 10; w = wrange(wi);
bc = betacrit(w,k,gamma,mu);

for ki = 1:length(krange)
    for wi = 1:length(wrange)
        if mod(wi,2)==0 && mod(ki,2)==0
            figure;
            positivity = valid(wi,ki,:);
            plot(brange,positivity(:))
            hold on;
            k = krange(ki); w = wrange(wi);
            bc = betacrit(w,k,gamma,mu);
            pl = plot([bc,bc],[-1,1]);
            pl.Color= 'black';
            pl.LineWidth = 2;
            ax = gca;
            ax.YLim = [min(positivity),max(positivity)];
            title(strcat('w =',num2str(wrange(wi)),' and k=',num2str(krange(ki))))    
        end
    end
end


% 
for bi = 1:length(brange)
    if mod(bi,2)==0
       b = brange(bi);
        figure;
        image(100*(valid(:,:,bi)>0))
        title(strcat('beta = ',num2str(b))); 
    end
end

bi = 15;
beta = brange(bi);

kcrits = @(beta,gamma,w,mu) (beta*(mu-1) - gamma - w);
kposs = ones(1,length(wrange));
kcrit = ones(1,length(wrange));
figure;
for wi = 1:length(wrange)
    valid_ws = valid(wi,:,bi);
    kposind = find(abs(valid_ws(:))<0.0001,1,'first');
    if isempty(kposind)
        kposs(wi) = NaN;
    else
        kposs(wi) = krange(kposind);
    end
    kcrit(wi) = kcrits(beta,gamma,wrange(wi),mu);
end
figure;
plot(krange,kposs)
hold on;
plot(krange,kcrit)



valid_ws = valid(wi,:,bi);
kposind = find(abs(valid_ws(:))<0.0001,1,'first');
if isempty(betacind)
    wposs(ki) = NaN;
else
    wposs(ki) = wrange(wposind);
end
plot(krange, valid_ws(:));



figure;
beta = 0.0025;
plot(krange, beta*(mu-1) - gamma - krange);
ax = gca; 
ax.XLim = [0,max(krange)];
ax.YLim = [0,max(krange)];
