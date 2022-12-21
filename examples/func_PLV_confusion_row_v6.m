function func_PLV_confusion_row_v6(PRCSDATADIR, DATASET, PERIOD, TIMING, FREQBAND)
% For given elec, calculate PLV values for all other electrodes in validE

%% Set Up variables and directories
% EEGLab won't call Octave functions in MATLAB, so load specific folders or
%       genpath and then rmpath the Octave functions
addpath('/home/despoB/hoycw/KnightLab/PRJ_DecisionMaking/Scripts/');
addpath(genpath('/home/knight/matar/MATLAB/toolbox/eeglab9_0_8_6b'));

PrcsDataDir = PRCSDATADIR;
uScorePos   = strfind(DATASET,'_'); periodPos = strfind(DATASET,'.');
SBJ         = DATASET(1:uScorePos-1);
task        = DATASET(uScorePos+1:periodPos-1);
curr_elec   = str2num(DATASET(periodPos+1:end));
WL          = 500;          % Window length for PLV in ms

%% Set Frequency and Timing parameters
% Load necessary variables
load(strcat(PrcsDataDir,'valid_elecs_allSBJandTasks.mat'));
validE = eval(['valid_elecs.' SBJ '.' task ';']);
SBJDataDir=strcat(PrcsDataDir,SBJ,'/',task,'/');
load(strcat(SBJDataDir,'gdat_notch.mat'));
load(strcat(SBJDataDir,'subj_globals.mat'),'srate');

switch FREQBAND
    case 'thetaL'
        loFreq = 3; hiFreq = 7;
    case 'theta'
        loFreq = 4; hiFreq = 8;
    case 'alpha'
        loFreq = 8; hiFreq= 12;
    case 'beta'
        loFreq = 15; hiFreq = 30;
    case {'cstULw','cstLw','cstHar','cstBeta'}
        [loFreq,hiFreq]=func_returnCstFBand(SBJ,task,FREQBAND);
    otherwise
        disp(strcat('Invalid frequency band variable: ',FREQBAND));
        error(strcat('Invalid frequency band variable: ',FREQBAND));
end

% Calculate windows of interest
WLinDP = WL*srate/1000;
[starts,ends]=func_returnWindows(SBJDataDir,task,PERIOD,TIMING,WLinDP);

% Prepare Directories of Results/Output
PLVvectorDir=strcat(SBJDataDir,'PLV_ePairs_trialsXtime/',FREQBAND,'_',PERIOD,TIMING,'/e',num2str(curr_elec,'%03d'),'/');
if ~isequal(exist(PLVvectorDir, 'dir'),7)
    mkdir(PLVvectorDir);
end

PLVconfRowDir = strcat(SBJDataDir,'PLV_confusionRows/',FREQBAND,'_',PERIOD,TIMING,'/');
PLVconfRowName=strcat(PLVconfRowDir,'PLV_',FREQBAND,'_',PERIOD,TIMING,'_e',num2str(curr_elec,'%03d'),'.mat');
if ~isequal(exist(PLVconfRowDir, 'dir'),7)
    mkdir(PLVconfRowDir);
end


%% PLV Calculations
% Allocate space for confusion matrix
PLV_electrodes = NaN(length(validE),1);
%PLV_pvalue = NaN(length(valid_elecs),1);

%calculate main elec data
elec1 = gdat(curr_elec,:);
elec1_highpass = eegfilt(elec1, srate, loFreq, []); %, 0, 3*fix(srate/loFreq), 0, 'fir1');
elec1_bandpass = eegfilt(elec1_highpass, srate, [], hiFreq); %, 0, 3*fix(srate/loFreq), 0, 'fir1');
hilbert1 = angle(hilbert(elec1_bandpass));

%calculate single PLV value and surrogate PLV for electrode pairs
for j = 1:length(validE)
    %calculate comparison elec data, including but not above the diagonal
    if (validE(j)<=curr_elec)
        elec2 = gdat(validE(j),:);
        elec2_highpass = eegfilt(elec2, srate, loFreq, []); %, 0, 3*fix(srate/loFreq), 0, 'fir1');
        elec2_bandpass = eegfilt(elec2_highpass, srate, [], hiFreq); %, 0, 3*fix(srate/loFreq), 0, 'fir1');
        hilbert2 = angle(hilbert(elec2_bandpass));
        
        %calculate PLV and pvalue for electrode pair
        PLV_electrodes(j) = func_PLV_singlePair(hilbert1, hilbert2, starts, ends, PLVvectorDir, curr_elec, validE(j));
        % PLV_pvalue(j) = PLV_surrogate(PLV_electrodes(j), hilbert1, hilbert2, start, finish, eventlength);
    else
        PLV_electrodes(j) = NaN;
    end
end

save(PLVconfRowName,'PLV_electrodes');

end
