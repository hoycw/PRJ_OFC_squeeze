function [surrogateBlobArea,surrogateBlobVolume,nativePLFvalues,...
    significantBlobAreaIndex,significantBlobVolumeIndex] = ...
    MonteCarloBlobCluster_PLF(data,trigger,window,nSurrs,...
    pBlobThresh,pRealThresh)
%[surrogateBlobArea,surrogateBlobVolume,nativePLFvalues,...
%     significantBlobAreaIndex,significantBlobVolumeIndex] = ...
%     MonteCarloBlobCluster_PLF(data,trigger,window,nSurrs,...
%     pBlobThresh,pRealThresh)
%
%MonteCarloBlobCluster_PLF ingests [time,frequency,channel] `data` matrix 
% of complex time-frequency transformed values with associated
% event `trigger` timestamps and generates a surrogate distribution of 2-D
% cluster statistics by circularly shifting these timestamps an `nSurrs` 
% number of times using a random (i.e. uniform) distribution of temporal
% offset jitters. The surrogate and original/native PLF time-frequency
% plots are thresholded using a within-frequency layer surrogate
% distribution and the quantile specified by `pBlobThresh`. Contiguous
% supra-threshold PLF values are grouped into clusters/blobs and a null 
% distribution of cluster statistics generated. Native PLF clusters are  
% then assessed using the `pRealThresh` level of significance.
%
%
%inputs:
% data: 2- or 3-D matrix of time-by-frequency-by-channel complex 
%   time-frequency values. data does NOT need to be normalized to have unit
%   length, the function will normalize all complex vector magnitudes to
%   have unit length.
% trigger: vector of integers indicating sample point indices corresponding
%   to event trigger times, range of [1,size(data,1)]
% window: vector of contiguous integers indicating sample point indices to 
%   include within the time analysis window, 0 indicates the exact time of 
%   the event trigger, <0 is prior to trigger, >0 is after trigger
% nSurrs: scalar positive integer indicating the number of random surrogate
%   shuffles to compute
% pBlobThresh: p-value used to enforce binary threshold on all surrogates 
%   and real data to create "blob" clusters, these blobs will then have 
%   their "area" (number of contiguous supra-threshold values) and "volume"
%   (summed PLF values within the contiguous cluster points) computed to 
%   generate a null distribution of surrogate cluster statistics, value 
%   should entered assuming a one-tailed test, range of (0,1)
% pRealThresh: p-value used to determine significance (i.e. alpha value) of
%   real data cluster statistics compared to distribution of surrogate 
%   cluster statistics, value should be entered under the assumpton it
%   is a one-tailed p-value, range of (0,1)
%
%outputs:
% surrogateBlobArea: cell array of dimensions [chan, 1] containing monte 
%   carlo generated distribution of blob cluster "area" (i.e. number of 
%   contiguous sample points)
% surrogateBlobVolume: cell array of dimensions [chan, 1] containing monte 
%   carlo generated distribution of blob cluster "volume" (i.e. sum of PLF 
%   values within each cluster)
% nativePLFvalues: PLF values of the originally input data-trigger
%   combination which can be subjected to thresholding and cluster 
%   statistic caluculations, of dimensions [numel(window), freq, chan]
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

%ensure window fits design constraints
if any(window~=round(window))
    warning(['values in `window` found to be ',...
        'not an integer. correction will be attempted.'])
    window=sort(unique(round(window)),'ascend');
end
window=reshape(window,[],1);

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
data=data(trigger(1)+window(1) : trigger(end)+window(end), :, :);
%update data variable size
ntp=size(data,1);
%shift time stamps to adjust for data trimming
trigger=trigger-(trigger(1)+window(1)-1);

%collect general size variables
nPulse=numel(trigger);
eventSamples=numel(window);

%create indexing variables
isEvent=false(ntp,1);
isEvent(bsxfun(@plus,window,trigger))=true;
%create surrogate time offsets
jitter=randi([1,ntp],[nSurrs,1]);

%initialize loop output variables
sorrogateBlobArea=cell(nch,1);
sorrogateBlobVolume=sorrogateBlobArea;
nativeAmpPvalues=nan(eventSamples,nfq,nch);
nativeAmpTvalues=nativeAmpPvalues;
significantAmpBlobAreaIndex=cell(nch,1);
significantAmpBlobVolumeIndex=significantAmpBlobAreaIndex;

surrogateBlobArea=cell(nch,1);
surrogateBlobVolume=surrogateBlobArea;
nativePLFvalues=nan(eventSamples,nfq,nch);
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
    chan=chan./abs(chan);
    %intialize parallel loop variables
    tempBlobArea=cell(nSurrs,1);
    tempBlobVolume=tempBlobArea;
    plfSurr=nan(eventSamples,nSurrs,nfq);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Monte Carlo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %parallel loop to creat PLF surrogates
    parfor nS=1:nSurrs
        %enact surrogate temporal shift
        shiftedChan=circshift(chan,jitter(nS),1);
        %compute surrogate phase-locking
        plfSurr(:,nS,:)=abs(mean(reshape(shiftedChan(isEvent,:),...
            eventSamples,nPulse,nfq),2));
    end
    %calculate surrogate frequency layer thresholds
    plfThresh=quantile(reshape(plfSurr,[],nfq),1-pBlobThresh);
    i=bsxfun(@ge,plfSurr,reshape(plfThresh,1,1,[]));
    %parallel loop through all surrogates to collect cluster stats
    parfor nS=1:nSurrs
        %slice surrogate
        plfNow=plfSurr(:,nS,:);
        %find contiguous significant p-values
        blob=bwconncomp(i(:,nS,:));
        %collect blob "surface area"
        tempBlobArea{nS}=cellfun(@(b) numel(b),blob.PixelIdxList);
        %collect blob "volumes"
        tempBlobVolume{nS}=cellfun(@(b) sum(plfNow(b)),...
            blob.PixelIdxList);
    end
    %store surrogate stats
    tempBlobArea=cat(2,tempBlobArea{:});
    surrogateBlobArea{nC}=tempBlobArea;
    tempBlobVolume=cat(2,tempBlobVolume{:});
    surrogateBlobVolume{nC}=tempBlobVolume;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%% calculate stat on original data %%%%%%%%%%%%%%%%
    %compute native phase-locking
    plf=squeeze(abs(mean(reshape(chan(isEvent,:),...
        eventSamples,nPulse,nfq),2)));
    nativePLFvalues(:,:,nC)=plf;
    %find contiguous significant values
    i=bsxfun(@ge,plf,plfThresh);
    blob=bwconncomp(i);
    %store significant blob area
    tempStat=cellfun(@(b) numel(b),blob.PixelIdxList);
    tempStat=blob.PixelIdxList(tempStat >= ...
        quantile(tempBlobArea,1-pRealThresh));
    significantBlobAreaIndex{nC}=unique(cat(1,tempStat{:}));
    %store significant blob area
    tempStat=cellfun(@(b) sum(plf(b)),blob.PixelIdxList);
    tempStat=blob.PixelIdxList(...
        tempStat>=quantile(tempBlobVolume,1-pRealThresh));
    significantBlobVolumeIndex{nC}=unique(cat(1,tempStat{:}));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %update user
    ['Monte Carlo simulation complete for ',num2str(nC),' of ',...
        num2str(nch),' channels after ',num2str(toc/60),'minutes']
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end