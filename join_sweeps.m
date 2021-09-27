%% join data sets

disp('loading datasets')
%% join the data sets from linux computer
% data set from linux

disp('loading sweepAprl17N500delta002.mat')
load('sweepAprl17N500delta002.mat')
mfaImax_17 = mfaImax;
mfaRinf_17 = mfaRinf;
simImax_17 = simImax;
simRinf_17 = simRinf;
% disp('the mfaRinf_17 at wi=1')
% disp(mfaRinf_17{1,1,1})
% disp('the mfaRinf_17 at wi=7')
% disp(mfaRinf_17{7,1,1})
% try
%     disp('the mfaRinf_17 at wi=8')
%     disp(mfaRinf_17{8,1,1})
% catch
%     disp('the mfaRinf_17 at wi=8 is not existing')
% end
% try
%     disp('the mfaRinf_17 at wi=9')
%     disp(mfaRinf_17{9,1,1})
% catch
%     disp('the mfaRinf_17 at wi=9 is not existing')
% end

% data set from windows with reverse wrange
disp('loading sweepAprl18N500delta002.mat')
load('sweepAprl18N500delta002.mat')
mfaImax_18 = mfaImax;
mfaRinf_18 = mfaRinf;
simImax_18 = simImax;
simRinf_18 = simRinf;
% disp('the mfaRinf_18 at wi=10')
% disp(mfaRinf_18{10,1,1})
% disp('the mfaRinf_18 at wi=11')
% disp(mfaRinf_18{11,1,1})
% try
%     disp(mfaRinf_18{22,1,1})
% catch
%     disp('at wi=22 there is no entry')
% end
    
 
% another data set from linux with correct wrange
disp('loading sweepAprl19N500delta002.mat')
load('sweepAprl19N500delta002.mat')
mfaImax_19 = mfaImax;
mfaRinf_19 = mfaRinf;
simImax_19 = simImax;
simRinf_19 = simRinf;
% disp('the mfaRinf_19 at wi=7')
% disp(mfaRinf_19{7,1,1})
% disp('the mfaRinf_19 at wi=8')
% disp(mfaRinf_19{8,1,1})
% disp('the mfaRinf_19 at wi=13')
% disp(mfaRinf_19{13,1,1})
% try
%     disp(mfaRinf_19{22,1,1})
% catch
%     disp('at wi=22 there is no entry')
% end

LW = length(wrange);

clear('mfaRinf')
clear('mfaImax')
clear('simRinf')
clear('simImax')

disp('start loading data.mat')
for wi=1:length(wrange)
    w=wrange(wi);
    for ki=1:length(kaprange)
        kap=kaprange(ki);
        for bi=1:length(brange)
            
            if wi<7
                mfaRinf{wi,ki,bi}=mfaRinf_17{wi,ki,bi};
                mfaImax{wi,ki,bi}=mfaImax_17{wi,ki,bi};
                simRinf{wi,ki,bi}=simRinf_17{wi,ki,bi};
                simImax{wi,ki,bi}=simImax_17{wi,ki,bi};
            elseif wi>=7 && wi<=13
                mfaRinf{wi,ki,bi}=mfaRinf_19{wi,ki,bi};
                mfaImax{wi,ki,bi}=mfaImax_19{wi,ki,bi};
                simRinf{wi,ki,bi}=simRinf_19{wi,ki,bi};
                simImax{wi,ki,bi}=simImax_19{wi,ki,bi};
            else
                mfaRinf{wi,ki,bi}=mfaRinf_18{wi,ki,bi};
                mfaImax{wi,ki,bi}=mfaImax_18{wi,ki,bi};
                simRinf{wi,ki,bi}=simRinf_18{wi,ki,bi};
                simImax{wi,ki,bi}=simImax_18{wi,ki,bi};
            end
            
            if bi==1
                disp(strcat('w = ',num2str(w),'. kap = ',num2str(kap),'.first entry of mfaRinf{wi,ki,1}', num2str(mfaRinf{wi,ki,bi})))
%                 disp(mfaRinf{wi,ki,bi})
            end
        end
    end
end

save(strcat('joinedsweep23AprN',num2str(N),'delta002.mat'),'mfaRinf','mfaImax','simRinf','simImax','wrange','kaprange','brange','N','mu','delta','Nt','tmax');
disp('done')