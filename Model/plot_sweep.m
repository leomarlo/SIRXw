% plot figures for each scenario

load('samplepaths22Jun.mat')
tzpes = {'infected','removed','si-links'};
SA =3;
same = 8;
Ni=1;
for N=Ns(2)
    ki=1;
    for kappa=kappas
        wi=1;
        for w=ws
            di=1;
            for delta=deltas(1)
                i0i=1;
                for i0=i0s
                    for sa=1:SA
                        ind_ind =1;
                        for ind=[3,4,6]
                            figure;
                            for sm=1:same
                                times= reshape(store_sm(Ni,ki,wi,i0i,sa,sm,1,:),Nt+1,[]);
                                values=reshape(store_sm(Ni,ki,wi,i0i,sa,sm,ind,:),Nt+1,[]);
                                plot(times,values/Ns(Ni),'LineWidth',1,'Color','blue') 
                                hold on;
                            end
                            times_mf = reshape(store_mf(Ni,ki,wi,i0i,sa,1,:),Nt+1,[]);
                            values_mf = reshape(store_mf(Ni,ki,wi,i0i,sa,ind,:),Nt+1,[]);
                            LT = find(times_mf>0,1,'last');
                            plot(times_mf(1:LT),values_mf(1:LT),'LineWidth',3,'Color','black')
                            titletxt = strcat(tzpes{ind_ind},' N=',num2str(N),' k=',num2str(kappa),' w=',num2str(w),' i0=',num2str(i0),' sa=',num2str(sa));
                            title(titletxt) 
                            resolution=300;
                            folder='figures/';
                            filename=strcat('Jun24_MF_vs_simulationN',num2str(N),'k',num2str(kappa),'w',num2str(w),'i0',num2str(i0),'sa',num2str(sa));
    %                         direction=strcat(folder,filename,'.png');
    %                         saveas(gcf,direction)
                            ind_ind = ind_ind+1;
                        end
                        
                       
                    end
                end
            end
        end
    end
end