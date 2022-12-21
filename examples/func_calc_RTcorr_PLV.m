function func_calc_RTcorr_PLV(PRCSDATADIR, DATASET, FREQBAND, sortID, pTHRESH)
%% Correlations between PLV on each trial and RTs
%   This version uses each window serparately

addpath('/home/despoB/hoycw/KnightLab/PRJ_DecisionMaking/Scripts/');
PrcsDataDir   = PRCSDATADIR;

auxName       = '';
pThresh       = str2num(pTHRESH);
TIMING        = {'Stim', 'Resp'};
PERIOD        = {'pre', 'post'};

uScorePos   = strfind(DATASET,'_');
SBJ         = DATASET(1:uScorePos-1);
task        = DATASET(uScorePos+1:end);
load(strcat(PrcsDataDir,'valid_elecs_allSBJandTasks.mat'));

%% Load SBJ information
SBJDataDir    = strcat(PrcsDataDir,SBJ,'/',task,'/');
ePairsDir     = strcat(SBJDataDir,'PLV_ePairs_trialsXtime/');
RTcorrDir     = strcat(SBJDataDir,'PLV_RTcorrelations/');
if ~isequal(exist(RTcorrDir, 'dir'),7)
    mkdir(RTcorrDir);
end
[startsName,~]=func_returnRTfiles(task);

% eInfo column contents:
% 1=elec#, 2=ROI, 3=Pattern, 4=Cluster, 5=idx in valid_elecs,
% 6=xLabel(e#), 7=yLabel(Patt-ROI), 8=Major divisions, 9=minor divisions
eInfoFileName = strcat(SBJDataDir,'eInfo_',sortID,'.mat');
load(eInfoFileName);
plot_elecs=[eInfo{:,1}];
load(strcat(SBJDataDir,'RTs_',startsName,'.mat'));

%% Calculations
numWins = length(TIMING)*length(PERIOD);
win=1;
corrMat = NaN(numWins,length(plot_elecs),length(plot_elecs));
pvalMat = NaN(numWins,length(plot_elecs),length(plot_elecs));
corrMat_thresh = NaN(numWins,length(plot_elecs),length(plot_elecs));
pvalMat_thresh = NaN(numWins,length(plot_elecs),length(plot_elecs));
for timing=1:length(TIMING)
    for period=1:length(PERIOD)
        % Load all PLV by Trails data, correlate with RT
        for newIdx1=1:length(plot_elecs)
            e1=eInfo{newIdx1,1};
            for newIdx2=1:length(plot_elecs)
                e2=eInfo{newIdx2,1};
                if e2<=e1
                    load(strcat(ePairsDir,FREQBAND,'_',PERIOD{period},TIMING{timing},auxName,...
                        '/e',num2str(e1,'%03d'),'/PLVe',num2str(e1,'%03d'),'_e',num2str(e2,'%03d'),'.mat'));
                    [corrMat(win,newIdx1,newIdx2),pvalMat(win,newIdx1,newIdx2)] = corr(PLV_vector,RTs);
                end
            end
        end
        corrMat(win,:,:)=func_fillOutMatrix(squeeze(corrMat(win,:,:)));
        pvalMat(win,:,:)=func_fillOutMatrix(squeeze(pvalMat(win,:,:)));
        
        % Mask by significance
        sigIdxs=find(pvalMat(win,:,:)<pThresh);
        tmp=zeros(size(squeeze(pvalMat(win,:,:))));
        tmp(sigIdxs)=1;
        pvalMat_thresh(win,:,:)=squeeze(pvalMat(win,:,:)).*tmp;
        corrMat_thresh(win,:,:)=squeeze(corrMat(win,:,:)).*tmp;
        
        win=win+1;
    end
end
corrName=strcat(SBJ,'_',task,'_RTcorr_',FREQBAND,'_combined',auxName,...
    '_',sortID,'_p',num2str(pThresh),'.mat');
save(strcat(RTcorrDir,corrName),'corrMat','pvalMat','corrMat_thresh','pvalMat_thresh');
end

