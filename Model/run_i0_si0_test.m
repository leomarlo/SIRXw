% run comparisin between initial infected and initial si links

N = 400;
mu = 15;
L = ceil( N * mu / 2 );

si0_range = (1/N):(1/N):(50/N);
i0_range = (1/N):(1/N):(5/N);
frequ = zeros(length(i0_range),length(si0_range)); 
samples = 10;
xsamps = 1000;
i0_ind = 1;
for i0 = i0_range
    for sa=1:samples
        A=ERG(N,L);
        for sx=1:xsamps
            x=InitState(i0,N);
            si0=((x==0)*A)*(x==1)'/N;
            si0_ind = find(si0==si0_range);
            frequ(i0_ind,si0_ind) =frequ(i0_ind,si0_ind)+1;
            frequency = num2str(int64(frequ(i0_ind,si0_ind)));
%             disp(frequency)
            switch frequency
                case '1'
                    save_As_1{i0_ind,si0_ind}=A;
                    save_xs_1{i0_ind,si0_ind}=x;
                case '2'
                    save_As_2{i0_ind,si0_ind}=A;
                    save_xs_2{i0_ind,si0_ind}=x;
                    disp(strcat('for i0=',num2str(i0),' and si0=',num2str(si0),' we already have ',num2str(frequ(i0_ind,si0_ind)),' entries.'))
                case '3'
                    save_As_3{i0_ind,si0_ind}=A;
                    save_xs_3{i0_ind,si0_ind}=x;
                case '4'
                    save_As_4{i0_ind,si0_ind}=A;
                    save_xs_4{i0_ind,si0_ind}=x;
                otherwise
                    u=0;
%                     disp(strcat('for i0=',num2str(i0),' and si0=',num2str(si0),' we already have ',num2str(frequ(i0_ind,si0_ind)),' entries.'))
            end
        end
        
    end
    i0_ind=i0_ind+1;
end
    
figure; image(frequ*10)

figure; plot(si0,frequ(2,:))

frequ
% 
lilsi.what='low_io_low_si0';
lilsi.si0=si0_range(14);
lilsi.i0=i0_range(2);
lilsi.A=save_As_1{2,14};
lilsi.x=save_xs_1{2,14};
lilsi.B=save_As_2{2,14};
lilsi.y=save_xs_2{2,14};
res{1}=lilsi;
hihsi.what='high_io_high_si0';
hihsi.si0=si0_range(48);
hihsi.i0=i0_range(5);
hihsi.A=save_As_1{5,48};
hihsi.x=save_xs_1{5,48};
hihsi.B=save_As_2{5,48};
hihsi.y=save_xs_2{5,48};
res{2}=hihsi;
lihsi.what='low_io_high_si0';
lihsi.si0=si0_range(48);
lihsi.i0=i0_range(2);
lihsi.A=save_As_1{2,48};
lihsi.x=save_xs_1{2,48};
lihsi.B=save_As_2{2,48};
lihsi.y=save_xs_2{2,48};
res{3}=lihsi;

% save('i0svssi0s.mat','res')
% 
% 

tmax=350;
Nt =1000;
beta=1/150;
gamma=1/40;
w=1/100;
rho0=1/100;
kap=12/10000;
delta=2/100;
pl =1;
SA =30;
for i=1:3
    A = res{i}.A;
    B = res{i}.B;
    x = res{i}.x;
    y = res{i}.y;
    clear('pl');
%     figure;

    Rinfs = zeros(1, 2*SA);
    for sa=1:SA
        disp('we are at sample number')
        disp(sa);
        resultA=SIRXi_w_fromAandx(A,x,beta,gamma,w,kap,delta,tmax,Nt); 
        Rinfs(2*(sa-1)+1)=max(resultA.NRs+resultA.NXis)/resultA.N;
%         pl(2*(sa-1)+1)=plot(resultA.times,((resultA.NRs+resultA.NXis)/resultA.N));
%         pl(2*(sa-1)+1).Color = 'black';
        resultB=SIRXi_w_fromAandx(B,y,beta,gamma,w,kap,delta,tmax,Nt);
        Rinfs(2*(sa-1)+2)=max(resultB.NRs+resultB.NXis)/resultB.N;
%         hold on;
%         pl(2*(sa-1)+2)=plot(resultB.times,((resultB.NRs+resultB.NXis)/resultB.N)); 
%         pl(2*(sa-1)+2).Color = 'blue';
%         disp(strcat('sa',num2str(sa),' i0=',num2str(res{i}.i0),' si0=',num2str(res{i}.si0),' sumx=',num2str(sum(res{i}.x))));
%         disp(strcat('sa',num2str(sa),' i0=',num2str(res{i}.i0),' si0=',num2str(res{i}.si0),' sumy=',num2str(sum(res{i}.y))));
    end
    res{i}.Rinfs = Rinfs;
    res{i}.Rinf_mean = mean(Rinfs);
    res{i}.Rinf_sdev = std(Rinfs);
%     ax =gca;
%     ax.XLabel.String ='Times';
%     ax.YLabel.String = 'rinfinity';
%     title(strcat('i0=',num2str(res{i}.i0),' and si0=',num2str(res{i}.si0)));
%     lgd = legend(pl);
    
end
    
%     result=SIRXi_w_fromAandx(Aini,x,beta,gamma,w,kap,delta,tmax,Nt);
%     
% end
for i=1:3
    disp([res{i}.i0,res{i}.si0,mean(res{i}.Rinfs(res{i}.Rinfs>=0.1))])
end