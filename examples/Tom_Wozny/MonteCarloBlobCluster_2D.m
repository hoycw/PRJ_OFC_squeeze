function [surrogateBlobArea,surrogateBlobVolume,...
    nativePvalues,nativeTvalues,...
    significantBlobAreaIndex,significantBlobVolumeIndex] = ...
    MonteCarloBlobCluster_2D(data,trigger,...
    prePulseWindow,postPulseWindow,nSurrs,pBlobThresh,pRealThresh)
%[surrogateBlobArea,surrogateBlobVolume,...
%     nativePvalues,nativeTvalues,...
%     significantBlobAreaIndex,significantBlobVolumeIndex] = ...
%     MonteCarloBlobCluster_2D(data,trigger,...
%     prePulseWindow,postPulseWindow,nSurrs,pBlobThresh,pRealThresh)
%
%MonteCarloBlobCluster_2D ingests [time,frequency,channel] `data` matrix 
% of time-frequency transformed amplitude values with associated
% event `trigger` timestamps and generates a surrogate distribution of 2-D
% cluster statistics by circularly shifting these timestamps an `nSurrs` 
% number of times using a random (i.e. uniform) distribution of temporal
% offset jitters. A T-test is computed across events, between the mean 
% `prePulseWindow` baseline and each 'postPulseWindow' timepoint. All
% surrogate and the native/original time-frequency plot are thresholded
% using the surrogate distribution of within-frequency layer T-values and a
% quantile threshold of `pBlobThresh`. Continguous supra-threshold values 
% are grouped into "blobs" or "clusters" and used to generate a null 
% distribution of blob areas (number of contiguous time points) and volumes
% (summed T-statistic within a blob), by which the native time triggered
% data threshold blobs can be compared using a surrogate alpha value of 
% `pRealThresh`.
%
%
%inputs:
% data: 2- or 3-D matrix of time-by-frequency-by-channel amplitude
%   time-frequency values. 
% trigger: vector of integers indicating sample point indices corresponding
%   to event trigger times, range of [1,size(data,1)]
% prePulseWindow: vector of integers indicating sample point 
%   indices to include for baseline mean comparisons prior to each event 
%   trigger, where 0 indicates the exact time of the event trigger, <0
%   prior to trigger, >0 after trigger
% postPulseWindow: vector of integers indicating sample point 
%   indices to include event related analysis, where 0 indicates the exact 
%   time of the event trigger, <0 prior to trigger, >0 after trigger
% nSurrs: scalar positive integer indicating the number of random surrogate
%   shuffles to compute
% pBlobThresh: p-value used to enforce binary threshold on all surrogates 
%   and real data to create "blob" clusters, these blobs will then have 
%   their "area" (number of contiguous supra-threshold values) and "volume"
%   (summed amplitude values within the contiguous cluster points) computed
%   to generate a null distribution of surrogate cluster statistics, value 
%   should entered as if a two-tailed test since increases and decreases in 
%   amplitude are being assessed, range of (0,1)
% pRealThresh: p-value used to determine significance (i.e. alpha value) of
%   real data cluster statistics compared to distribution of surrogate 
%   cluster statistics, value should entered as if a one-tailed test 
%   (one-tailed p-value will be used to compare blob areas since only 
%   looking for larger areas, pRealThresh will automatically be cut in half
%   to be used as a two-tailed value for volume comparisons since increases
%   and decreases are being assessed), range of (0,1)
%
%outputs:
% surrogateBlobArea: cell array of dimensions [chan, 1] containing monte 
%   carlo generated distribution of blob cluster "area" (i.e. number of 
%   contiguous sample points)
% surrogateBlobVolume: cell array of dimensions [chan, 1] containing monte 
%   carlo generated distribution of blob cluster "volume" (i.e. sum of PLF 
%   values within each cluster)
% nativePvalues: T-test p-values of the originally input data-trigger
%   combination which can be subjected to thresholding and cluster 
%   statistic caluculations, size [numel(postPulseWindow), freq, chan]
% nativeTvalues: T-statistic values of the originally input data-trigger
%   combination which can be sujected to thresholding and cluster statistic
%   caluculations, size [numel(postPulseWindow), freq, chan]
% significantBlobAreaIndex: cell array of dimensions [chan, 1] containing 
%   per-channel indices of significant blob area sample points within the
%   `window` time frame provided
% significantBlobVolumeIndex: cell array of dimensions [chan, 1] containing 
%   per-channel indices of significant blob volume sample points within the
%   `window` time frame provided
%
%
% by Thomas Wozny
% last update 12/10/22


%%%%%%%%%%%%%%%%%%%%%%%% initialize variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data variable size
[ntp,nfq,nch]=size(data);

%ensure trigger fit design constraints
if any(trigger<1) || any(trigger>ntp) || any(trigger~=round(trigger))
    warning(['values in `trigger` found to be out of the range ',...
        '[1, size(data,1)] and/or not an integer. '...
        'correction will be attempted.'])
    trigger=sort(unique(abs(round(trigger))),'ascend');
end
trigger=reshape(trigger,1,[]);

%ensure peri-trigger windows fit design constraints
if any(prePulseWindow~=round(prePulseWindow))
    warning(['values in `prePulseWindow` found to be ',...
        'not an integer. correction will be attempted.'])
    prePulseWindow=sort(unique(round(prePulseWindow)),'ascend');
end
prePulseWindow=reshape(prePulseWindow,[],1);
if any(postPulseWindow~=round(postPulseWindow))
    warning(['values in `postPulseWindow` found to be ',...
        'not an integer. correction will be attempted.'])
    postPulseWindow=sort(unique(round(postPulseWindow)),'ascend');
end
postPulseWindow=reshape(postPulseWindow,[],1);

%ensure # of surrogates fit design constraints
if numel(nSurrs)~=1 || any(nSurrs~=round(nSurrs)) || any(nSurrs<=0)
    warning(['values in `nSurrs` found to be </=0, not an integer, ',...
        'and/or not a scalar. correction will be attempted.'])
    nSurrs=abs(round(max(nSurrs)));
end

%ensure p-values fit design constraints
if numel(pBlobThresh)~=1 ||  any(pBlobThresh<=0) || any(pBlobThresh>=1)
    warning(['values in `pBlobThresh` found to be out of the range ',...
        '(0,1) and/or not a scalar. correction will be attempted.'])
    pBlobThresh=abs(pBlobThresh);
    pBlobThresh=min(pBlobThresh(pBlobThresh>0 & pBlobThresh<1));
end
if numel(pRealThresh)~=1 ||  any(pRealThresh<=0) || any(pRealThresh>=1)
    warning(['values in `pRealThresh` found to be out of the range ',...
        '(0,1) and/or not a scalar. correction will be attempted.'])
    pRealThresh=abs(pRealThresh);
    pRealThresh=min(pRealThresh(pRealThresh>0 & pRealThresh<1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%% variable housekeeping prior to loop %%%%%%%%%%%%%%%%%%%%
%trim data from excessive time points
data=data(trigger(1)+prePulseWindow(1) : ...
    trigger(end)+postPulseWindow(end), :, :);
%update data variable size
ntp=size(data,1);
%shift time stamps to adjust for data trimming
trigger=trigger-(trigger(1)+prePulseWindow(1)-1);

%collect general size variables
nPulse=numel(trigger);
baseSamples=numel(prePulseWindow);
eventSamples=numel(postPulseWindow);

%create indexing variables
isBase=false(ntp,1);
isEvent=isBase;
isBase(bsxfun(@plus,prePulseWindow,trigger))=true;
isEvent(bsxfun(@plus,postPulseWindow,trigger))=true;
%create surrogate time offsets
jitter=randi([1,ntp],[nSurrs,1]);

%initialize loop output variables
sorrogateBlobArea=cell(nch,1);
sorrogateBlobVolume=sorrogateBlobArea;
nativePvalues=nan(eventSamples,nfq,nch);
nativeTvalues=nativePvalues;
significantBlobAreaIndex=cell(nch,1);
significantBlobVolumeIndex=significantBlobAreaIndex;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%% master channel loop %%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop through each channel
for nC=1:nch
    %start channel timer
    tic;
    %grab current channel and normalize
    chan=data(:,:,nC);
    %intialize parallel loop variables
    tempBlobArea=cell(nSurrs,1);
    tempBlobVolume=tempBlobArea;
    ampTsurr=nan(eventSamples,nSurrs,nfq);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Monte Carlo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %parallel loop to create amplitude T-stat surrogates
    parfor nS=1:nSurrs
        %enact surrogate temporal shift
        shiftedChan=circshift(chan,jitter(nS),1);
        %gather baseline for comparison
        base=reshape(mean(reshape(shiftedChan(isBase,:),...
            baseSamples,nPulse,nfq),1),nPulse,1,nfq);
        %T test
        [~,~,~,stat]=ttest(repmat(base,1,eventSamples,1),permute(reshape...
            (shiftedChan(isEvent,:),eventSamples,nPulse,nfq),[2,1,3]));
        ampTsurr(:,nS,:)=squeeze(stat.tstat);
    end
    %calculate surrogate frequency layer thresholds
    ampThresh=quantile(reshape(ampTsurr,[],nfq),...
        [pBlobThresh,1-pBlobThresh]);
    i=bsxfun(@le,ampTsurr,reshape(ampThresh(1,:),1,1,[])) | ...
        bsxfun(@ge,ampTsurr,reshape(ampThresh(2,:),1,1,[]));
    %parallel loop through all surrogates to collect cluster stats
    parfor nS=1:nSurrs
        %slice surrogate
        ampNow=ampTsurr(:,nS,:);
        %find contiguous significant p-values
        blob=bwconncomp(i(:,nS,:));
        %collect blob "surface area"
        tempBlobArea{nS}=cellfun(@(b) numel(b),blob.PixelIdxList);
        %collect blob "volumes"
        tempBlobVolume{nS}=cellfun(@(b) sum(ampNow(b)),...
            blob.PixelIdxList);
    end
    %store surrogate stats
    tempBlobArea=cat(2,tempBlobArea{:});
    surrogateBlobArea{nC}=tempBlobArea;
    tempBlobVolume=cat(2,tempBlobVolume{:});
    surrogateBlobVolume{nC}=tempBlobVolume;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%% calculate stat on original data %%%%%%%%%%%%%%%%
    %grab baseline mean
    base=reshape(mean(reshape(chan(isBase,:),...
        baseSamples,nPulse,nfq),1),nPulse,1,nfq);
    %T test
    [~,p,~,stat]=ttest(repmat(base,1,eventSamples,1), permute(...
        reshape(chan(isEvent,:),eventSamples,nPulse,nfq),[2,1,3]));
    p=squeeze(p);
    T=squeeze(stat.tstat);
    %stor native stats
    nativePvalues(:,:,nC)=p;
    nativeTvalues(:,:,nC)=T;
    %find contiguous significant p-values
    blob=bwconncomp(bsxfun(@le,T,ampThresh(1,:)) | ...
        bsxfun(@ge,T,ampThresh(2,:)));
    %store significant blob area
    tempStat=cellfun(@(b) numel(b),blob.PixelIdxList);
    tempStat=blob.PixelIdxList(tempStat >= ...
        quantile(tempBlobArea,1-pRealThresh));
    significantBlobAreaIndex{nC}=unique(cat(1,tempStat{:}));
    %store significant blob area
    tempStat=cellfun(@(b) sum(T(b)),blob.PixelIdxList);
    tempStat=blob.PixelIdxList(...
        tempStat<=quantile(tempBlobVolume,pRealThresh*0.5) | ...
        tempStat>=quantile(tempBlobVolume,1-pRealThresh*0.5));
    significantBlobVolumeIndex{nC}=unique(cat(1,tempStat{:}));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %update user
    ['Monte Carlo simulation complete for ',num2str(nC),' of ',...
        num2str(nch),' channels after ',num2str(toc/60),'minutes']
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end