function [M,info,I,Q,clipping]=Decode_Acqdecode_Raw(subj,dateA,suf,params)
%
% This function loads and decodes raw data from Acqdecode
% The function must be run in the folder containing the data files.
% The function first reads the Acqdecode-out info file.
% Then the function reads the *.cfg file for the encoding pattern. This
% file should be in the .../Support_files/ folder, possibly under a 
% .../cfg files/ folder.
% params is a structure that can contain the following parameters:
%   decode_type     {'IQ','FFT'} decoding strategies. 'FFT' is default
%   fftNh           if 'FFT' is decode type, fftNh selects how many
%                   harmonics are used in decoding. (1 or 2)
%   nfft            {'np2','spts'} Selects mode for the number of time points 
%                   to use in the fft of each time step. 'spts' is default
%   Wtype           window type. 'Hann' is default
%   adcT.i          ADC transient kill start point. default=10
%   adcT.f          ADS transient kill end point. default=4
%   noMatFile       Dsiable saving .mat files at end
%   fftIntPeriods   Sets FFT length to approximately an integer number of
%                   periods for each frequency
% Infrastructure for all systems utilizes the maximum allowed A/D channel
% per computer is the frame synch and the 2nd-to-last channel is the stim
% synch.



%% Parameters and Initialization
mPHeight=0.1;
mPDist=1000;
ext='.raw';
if ~exist('params','var'),params=[];end
if ~isfield(params,'overwrite'),params.overwrite=1;end
if ~isfield(params,'sr'),params.sr=96e3;end
sr=params.sr;
if ~isfield(params,'ND2'),params.ND2=0;end
if ~isfield(params,'tag')
    datafile2save=[dateA,'-',subj,'-',suf]; % Data file name
else
    datafile2save=[dateA,'-',subj,'-',suf,'_',params.tag]; % Data file name
end
datafile=[dateA,'-',subj,'-',suf]; % Data file name

if exist([datafile2save,'.mat'],'file') && ~params.overwrite
    load([datafile2save,'.mat'])
    if isfield(info, 'MEAS')
        if istable(info.MEAS)
            info.MEAS = table2struct(info.MEAS, 'ToScalar', true)
        end
    end
    if isfield(info, 'pairs')
        if istable(info.pairs)
            info.pairs = table2struct(info.pairs, 'ToScalar', true)
        end
    end
else % do all of decoding...
    
if ~isfield(params,'decode_type') % options: 'IQ','FFT'
    params.decode_type='FFT'; % Changed default to FFT 180131
end 
if ~isfield(params,'fftNh'),params.fftNh=1;end % options: 1,2
if ~isfield(params,'nfft') % options: 'np2','spts'
    params.nfft='spts'; % Changed default to spts 180213
end 
if ~isfield(params,'Wtype') % 
    params.Wtype='Hann'; %  
end 
if ~isfield(params,'adcT') % Enforce adc transient 180213
    params.adcT.i=10; % initial dead window (10?) (12?)
    params.adcT.f=4; % ending dead window (4?) (5?)
    t0=params.adcT.i;
else
    t0=0;
end 
if ~isfield(params,'clipM') % clipping max
    params.clipM=1; 
    params.clipMu=0.5; 
end 

% Also, you can have +/- thresholds if necessary
if ~isfield(params,'clipPlus') % clippling+
    params.clipPlus = +params.clipM ;
end
if ~isfield(params,'clipMinus')% clippling-
    params.clipMinus = -params.clipM ;
end
if ~isfield(params,'clipRange')%clipping range
    params.clipRange = 1 ;
end
if ~isfield(params,'New2pass')%clipping range
    params.New2pass = 1 ;
end

if ~isfield(params,'noMatFile') % disable mat file saving
    params.noMatFile = 0 ;
end
if ~isfield(params,'fftIntPeriods') % enforce integer fft periods
    params.fftIntPeriods = 0 ;
end


%% Load info file
info=readkeyfile([datafile,'-info.txt']);
Ns=info.srcnum;
Nd=info.detnum;
Nwl=info.cnum;
Nts=info.nts;
switch length(num2str(info.run))
    case 1
fn_base=[dateA,'-run00',num2str(info.run),'-ch'];
    case 2
fn_base=[dateA,'-run0',num2str(info.run),'-ch'];
end


if isfield(info, 'pairs')
    if istable(info.pairs)
        info.pairs = table2struct(info.pairs, 'ToScalar', true)
    end
end
if ~isfield(params,'FSchan')
switch info.pad % These must show the max A/D channel for each system
    case 'Adult_96x92a'
        FSchan=48; % info.nmotu * #chan in Motu (12)
    case 'Adult_96x92b'
        FSchan=48; % info.nmotu * #chan in Motu (12)
    case 'AdultV24x28'
        FSchan=28; % info.nmotu * #chan in Motu (12)
    case 'Baby32x34'
        FSchan=36; % info.nmotu * #chan in Motu (12)
    case 'AdultClinical2_48x34'
        FSchan=48; % info.nmotu * #chan in RME (32)
    case 'Matador_84x91'
        FSchan=96; % info.nmotu * #chan in FocusRite (16)
    case 'UHD126x126'
        FSchan=128; % Max focusrite channels as of 180705
    case 'Gates_96x92'
        FSchan=128; % Max focusrite channels as of 180705
    case 'Gates_128x125'
        FSchan=128; % Max focusrite channels as of 180705
    case 'Matador_255x252_201009_half_a'
        FSchan=128; % Max focusrite channels as of 180705
    case 'Matador_255x252_201009_half_b'
        FSchan=128; % Max focusrite channels as of 180705
    case 'Adult_128x125'
        FSchan=128; % Max focusrite channels as of 180705
end
else
    FSchan=params.FSchan;
end


%% Load encoding pattern cfg file
EncCfg=readkeyfile([info.enc,'.cfg']);
Nreg=EncCfg.nreg;
NregReps=ceil(info.srcnum./(info.nts*Nreg));
info.EncCfg=EncCfg;
flashPerc=EncCfg.nisamp/(EncCfg.nisamp+EncCfg.niblank);
periodSamp=sr./EncCfg.freq; % Number of samples per period for ea freq
Nfreq=length(EncCfg.freq);

if ~isfield(params,'FSchan')
    if isfield(EncCfg,'adcmax'), FSchan=EncCfg.adcmax;end
end

SSchan=FSchan-1; % stim synch (may not be used for now)

%% Get frame synch raw file and Calculate empirical sampling info
disp('<<< Finding frame synch points')
[rawFilename,rawFilenameError] = fnbase2rawFileName( fn_base, FSchan, ext ) ;
fid=fopen(rawFilename) ;

Fsynch=fread(fid,'float32');fclose(fid);

% Find frame synchs
if abs(min(Fsynch))>abs(max(Fsynch)) % Adjust for possible neg synch pulse
    Fsynch=-Fsynch;
end

DtFsynch=-([Fsynch(1);Fsynch]-[Fsynch;Fsynch(end)]);
DtFsynch(end)=[];
[~,fs]=findpeaks(DtFsynch,'MinPeakHeight',mPHeight,'MinPeakDistance',mPDist);
Nf=length(fs)-1;
clear DtFsynch Fsynch

% Establish ON sampling time for each frequency
if length(unique(diff(fs)))>2
    disp(['Warning, variable frame sizes; ',num2str(unique(diff(fs(:)))'),...
        'frames may have been dropped'])
    figure;histogram((diff(fs)),[min(diff(fs))-1:max(diff(fs))+1]);
    set(gca,'YScale','log');
    xlabel('Samples per frame');ylabel('Number frames')
    pause(2)
%     return
end
% Flength=round(mean(diff(fs)));  % Assuming no frame synchs were missed
Flength=mode(diff(fs));  % Assuming no frame synchs were missed
disp(['<<< Frame rate ~ ',num2str(sr/Flength),' Hz'])

info.nframe=length(fs)-1; % number of full acquired frames
info.framesize=Flength; % mean frame length
info.framerate=sr/Flength; % calculate empirical frame rate (MOTU @ 96kHz)

SamplesPerTimeStepMax=floor((Flength/(EncCfg.nts))*flashPerc)-...
    params.adcT.i-params.adcT.f ;
numPeriods=floor(SamplesPerTimeStepMax./periodSamp);
SamplesPerTimeStepIntPeriods=floor(numPeriods.*periodSamp); 
if ( params.fftIntPeriods == 1 )
    SamplesPerTimeStep = SamplesPerTimeStepIntPeriods ;
else
    SamplesPerTimeStep = ones(size(EncCfg.freq)).*SamplesPerTimeStepMax ;
end

% prep for fft-based decoding
switch params.nfft
    case 'np2' % Changed default to FFT 180131
        Ndft = pow2(nextpow2(SamplesPerTimeStep));   
    case 'spts'
    Ndft = SamplesPerTimeStep;   % DO NOT Zero pack to a power of 2. (ER)
end
for j=1:Nfreq % indices associated w mod freqs
    freqs=sr*(0:(Ndft(j)/2))./Ndft(j);    % domain of FFT
    [mF,idxf]=min(abs(freqs-EncCfg.freq(j))); % Find FFT bin closest to EncCfg.freq(i)
    switch params.fftNh
        case 1
            idxF{j}=[idxf-1,idxf,idxf+1];
        case 2
            [mF,idxfh]=min(abs(freqs-2*EncCfg.freq(j)));
            idxF{j}=[idxf-1,idxf,idxf+1,idxfh-1,idxfh,idxfh+1];
    end
end

params.SamplesPerTimeStep=SamplesPerTimeStep;
info.DecodingParams=params;

%% Set up Hamm window and sin/cos basic functions for ea freq
[HammF,sinF,cosF]=Generate_Decoding_Bases(...
    SamplesPerTimeStep,sr,EncCfg.freq,params);


%% Decode from raw data
disp('<<< Decoding raw data')
M=zeros(Ns,Nd,Nwl,Nf);
clipping=zeros(Ns,Nd,Nf);
if strcmp(params.decode_type,'IQ')
I=zeros(Ns,Nd,Nwl,Nf);
Q=zeros(Ns,Nd,Nwl,Nf);
else
   I=[];
   Q=[]; 
end

for d=1:Nd
    % Load data
    fprintf('d%s ',num2str(d))
    [rawFilename,rawFilenameError]=fnbase2rawFileName( fn_base, d, ext );
    fid=fopen(rawFilename);
    APDdata=fread(fid,'float32');fclose(fid);
    
    % Choose strategy here
    switch params.decode_type
        case 'IQ'
    
    % Decode all frames at once: IQ
    for ts=1:Nts 
        % For a given time step, find start times for ea frame
        ti=round(fs+(ts-1)*diff([fs;fs(end)])./Nts);ti(end)=[];
        
        tIdx=bsxfun(@plus,ti,[t0:(t0+SamplesPerTimeStep(1)-1)])';
        d1=reshape(APDdata(tIdx(:)),SamplesPerTimeStep(1),[])';
        clipped=DetectClipping(d1,params)';
        for n=1:Nreg
            if ( ts+((n-1)*Nts) <= info.srcnum )
                clipping(ts+((n-1)*Nts),d,:)=clipped;
            end
        end
            
        [M(ts,d,1,:),I(ts,d,1,:),Q(ts,d,1,:)]=MIQ_Decode_Raw_Acqdecode(...
            d1,sinF{1},cosF{1},HammF{1},EncCfg.div(1));% REGION 1, WL 1
        
        tIdx=bsxfun(@plus,ti,[t0:(t0+SamplesPerTimeStep(2)-1)])';
        d1=reshape(APDdata(tIdx(:)),SamplesPerTimeStep(2),[])';
        [M(ts,d,2,:),I(ts,d,2,:),Q(ts,d,2,:)]=MIQ_Decode_Raw_Acqdecode(...
            d1,sinF{2},cosF{2},HammF{2},EncCfg.div(2));% REGION 1, WL 2
        
        if Nreg>1
        tIdx=bsxfun(@plus,ti,[t0:(t0+SamplesPerTimeStep(3)-1)])';
        d1=reshape(APDdata(tIdx(:)),SamplesPerTimeStep(3),[])';
        [M(ts+Nts,d,1,:),I(ts+Nts,d,1,:),Q(ts+Nts,d,1,:)]=...
            MIQ_Decode_Raw_Acqdecode(...
            d1,sinF{3},cosF{3},HammF{3},EncCfg.div(3));% REGION 2, WL 1
        
        tIdx=bsxfun(@plus,ti,[t0:(t0+SamplesPerTimeStep(4)-1)])';
        d1=reshape(APDdata(tIdx(:)),SamplesPerTimeStep(4),[])';
        [M(ts+Nts,d,2,:),I(ts+Nts,d,2,:),Q(ts+Nts,d,2,:)]=...
            MIQ_Decode_Raw_Acqdecode_U87(...
            d1,sinF{4},cosF{4},HammF{4},EncCfg.div(4));% REGION 2, WL 2
        end
    end
    
    
        case 'FFT'
    % Decode all frames at once: Full FFT
    for ts=1:Nts 
        % For a given time step, find start times for ea frame
        ti=round(fs+(ts-1)*diff([fs;fs(end)])./Nts);ti(end)=[];        

        if ( params.fftIntPeriods == 1 )
            % Number of samples in FFT different for each frequency
            tIdx=bsxfun(@plus,ti,[t0:(t0+SamplesPerTimeStep(1)-1)])';
            d1=reshape(APDdata(tIdx(:)),SamplesPerTimeStep(1),[])';
            clipped=DetectClipping(d1,params)'; 
            for n=1:Nreg
                if ( ts+((n-1)*Nts) <= info.srcnum )
                    clipping(ts+((n-1)*Nts),d,:)=clipped;
                end
            end
            
            M(ts,d,1,:)=FFT_Decode_Raw_Acqdecode(...
                d1,HammF{1},Ndft(1),idxF{1},EncCfg.div(1));

            tIdx=bsxfun(@plus,ti,[t0:(t0+SamplesPerTimeStep(2)-1)])';
            d1=reshape(APDdata(tIdx(:)),SamplesPerTimeStep(2),[])';
            M(ts,d,2,:)=FFT_Decode_Raw_Acqdecode(d1,HammF{2},...
                Ndft(2),idxF{2},EncCfg.div(2));

            if Nreg>1
                tIdx=bsxfun(@plus,ti,[t0:(t0+SamplesPerTimeStep(3)-1)])';
                d1=reshape(APDdata(tIdx(:)),SamplesPerTimeStep(3),[])';
                M(ts+Nts,d,1,:)=FFT_Decode_Raw_Acqdecode(...
                    d1,HammF{3},Ndft(3),idxF{3},EncCfg.div(3));

                tIdx=bsxfun(@plus,ti,[t0:(t0+SamplesPerTimeStep(4)-1)])';
                d1=reshape(APDdata(tIdx(:)),SamplesPerTimeStep(4),[])';
                M(ts+Nts,d,2,:)=FFT_Decode_Raw_Acqdecode(...
                    d1,HammF{4},Ndft(4),idxF{4},EncCfg.div(4));
            end
        else
            % Number of samples in FFT same for each frequency
            tIdx=bsxfun(@plus,ti,[t0:(t0+SamplesPerTimeStep(1)-1)])';
            d1=reshape(APDdata(tIdx(:)),SamplesPerTimeStep(1),[])';
%             clipping(ts:Nts:end,d,:)=repmat(DetectClipping(d1,params)',Nreg,1,1); 
            clipped=DetectClipping(d1,params)'; 
            mags = FFT_Decode_Raw_Acqdecode(...
                d1,HammF{1},Ndft(1),idxF,EncCfg.div);
            
            for n=1:Nreg
                m0=((n-1)*Nwl)+1;
                mF=n*Nwl;
                for Nrr=1:NregReps % better repeated region support 190325
                    sources=ts+((n-1)*Nts)+(Nrr-1)*(Nreg*Nts);
                    if ( sources <= info.srcnum )
                        M(sources,d,:,:) = mags(m0:mF,:);
                        clipping(sources,d,:)=clipped;
                    end
                end
            end
        end       
    end
    
    end
end
clear APDdata


%% If 2-pass, fix data and frames before synch
% Idea is to keep hot data overall and replace clipped hot data with cool
% break into even and odd
% find hot half per region
% make full hot set
% find clipped in hot and replace with cool set
if strfind(info.enc,'2pass')
    
    if params.New2pass
        Nf_new=floor(size(M,4)/2); % Number of frames after merging hot and cool        
        oM=M(:,:,:,1:2:(Nf_new*2));
        eM=M(:,:,:,2:2:(Nf_new*2));
        dMmu=mean(mean(eM-oM,4),3);
        Mhot=zeros(size(oM));
        Mcool=zeros(size(oM));
        clippingHot=zeros(Ns,Nd,Nf_new);
        clippingCool=zeros(Ns,Nd,Nf_new);
        OddRegionHotness=zeros(Nreg,1);
        for j=1:Nreg
            si=Nts*(j-1)+1;
            sf=Nts*j;
            if sf>info.srcnum, sf=info.srcnum;end
            clippingRegion = clipping(si:sf,:,(1:Nf_new*2));
        
            %%%% Find Region-specific hot and cool frames UPDATED 20190311
            OddRegionHotness(j)=sum(sum(dMmu(si:sf,:),2))<0; 
            if j>1
                if strfind(info.enc,'Gates')
                    OddRegionHotness(j)=sum(sum(dMmu(si:sf,:),2))<0;                    
                else
                    OddRegionHotness(j)=~OddRegionHotness(j-1);                    
                end
%                 
%                 if EncCfg.freq(1)==EncCfg.freq(1+EncCfg.cnum)
%                     OddRegionHotness(j)=~OddRegionHotness(j-1);
%                 else
%                     OddRegionHotness(j)=sum(sum(dMmu(si:sf,:),2))<0;
%                 end
            end
            
            if OddRegionHotness(j)
                dHot=oM(si:sf,:,:,:);
                dCool=eM(si:sf,:,:,:);
                clippingRegionHot=clippingRegion(:,:,1:2:end);
                clippingRegionCool=clippingRegion(:,:,2:2:end);
                fsKeep=fs(1:2:end);
                disp(['<<<Interpolating odd/even bright/dim data',...
                    ' for Region ',num2str(j)])
                info.system.framerate=info.framerate;
                dHrs=resample_tts(dHot,info,2*info.framerate);
                dHot=dHrs(:,:,:,2:2:end);
                dHot(dHot<0)=eps;
                clear dHrs
            else
                dHot=eM(si:sf,:,:,:);
                dCool=oM(si:sf,:,:,:);
                clippingRegionHot=clippingRegion(:,:,2:2:end);
                clippingRegionCool=clippingRegion(:,:,1:2:end); 
                fsKeep=fs(2:2:end);
                disp(['<<<Interpolating even/odd bright/dim data',...
                    ' for Region ',num2str(j)])
                info.system.framerate=info.framerate;
                dCrs=resample_tts(dCool,info,2*info.framerate);
                dCool=dCrs(:,:,:,2:2:end);
                dCool(dCool<0)=eps;
                clear dCrs                
            end
            
            %%%% Get Scaling for dc cool to hot and scale cool
            Nset=(Nreg*Nwl);
            DimScale=zeros(Nwl,1);
            for wl=1:Nwl
                DimScale(wl)=[sin(pi*info.EncCfg.dc(wl+(j-1)*Nwl))]/...
                    [sin(pi*info.EncCfg.dc(wl+Nset+(j-1)*Nwl))];
                DimScale(wl)=max(DimScale(wl),1/DimScale(wl));
                dCool(:,:,wl,:)=dCool(:,:,wl,:).*DimScale(wl);
            end
            
            % Re-populate temporary arrays
            Mhot(si:sf,:,:,:)=dHot;
            Mcool(si:sf,:,:,:)=dCool;
            clippingHot(si:sf,:,:)=clippingRegionHot; 
            clippingCool(si:sf,:,:)=clippingRegionCool;            
        end
        
        % Put Cool data in with hot data
        Mhot=reshape(Mhot,Ns*Nd,Nwl,[]);
        Mcool=reshape(Mcool,Ns*Nd,Nwl,[]);
        clippingHot=reshape(clippingHot,Ns*Nd,[]);
        clippingCool=reshape(clippingCool,Ns*Nd,[]);
        
        Cool2Hot=sum(clippingHot,2)>0;
        Mhot(Cool2Hot,:,:)=Mcool(Cool2Hot,:,:);
        clippingHot(Cool2Hot,:)=clippingCool(Cool2Hot,:);
        
        MBestPassAllRegions = reshape(Mhot,Ns,Nd,Nwl,[]);
        clippingAllRegions = reshape(clippingHot,Ns,Nd,[]);
        
    
        
        
    else
    
    for j = 0:Nreg-1
        % Get M, clipping for each region
        RegSrcStartIdx = j*Nts + 1 ;
        RegSrcEndIdx = (j+1)*Nts ; 
        if ( RegSrcEndIdx > info.srcnum )
            RegSrcEndIdx = info.srcnum ;
        end
        Nf_new=floor(size(M,4)/2); % Number of frames after merging hot and cool
        % Throw out last frame if there is an odd number so that 
        % dHot and dCool are the same size
        MRegion = M(RegSrcStartIdx:RegSrcEndIdx,:,:,1:(Nf_new*2)); 
        clippingRegion = clipping(RegSrcStartIdx:RegSrcEndIdx,:,(1:Nf_new*2)) ;
        
        % Scale each wavelengths data to correct for dc val
        % scaling=[sin(pi * high dc)]/[sin(pi * low dc)]
%       % indices must match Region/freq dc vals
%         DimScale(1)=sin(pi*info.EncCfg.dc(1+4*j))/sin(pi*info.EncCfg.dc(3+4*j));
%         DimScale(2)=sin(pi*info.EncCfg.dc(2+4*j))/sin(pi*info.EncCfg.dc(4+4*j));
        % We want the large ratio because we are going to scale up the
        % lower duty cycle signal but we don't know which is which yet        
%         DimScale(1)=max(DimScale(1),1/DimScale(1)) ;
%         DimScale(2)=max(DimScale(2),1/DimScale(2)) ;
        
        Nset=(Nreg*Nwl);
        DimScale=zeros(Nwl,1);
        for wl=1:Nwl
            DimScale(wl)=[sin(pi*info.EncCfg.dc(wl+j*Nwl))]/...
                [sin(pi*info.EncCfg.dc(wl+Nset+j*Nwl))];
            DimScale(wl)=max(DimScale(wl),1/DimScale(wl));
        end
            

        
        % Only look for clipping for S/D dist < 4 cm to determine which 
        % pass is clipped because it could be
        % caused by bright sources in another region
        
        
        Rad=load(['Pad_',info.pad]); % Load Pad file for cap
        radR2dSrcByDet = reshape(Rad.info.pairs.r2d(1:(Ns*Nd)),Ns,Nd);
        
        
        radR2dRegionSrcByDet = radR2dSrcByDet(RegSrcStartIdx:RegSrcEndIdx,:) ;
        radR2dRegionSrcXDet = reshape(radR2dRegionSrcByDet,...
            size(radR2dRegionSrcByDet,1)*size(radR2dRegionSrcByDet,2),1) ;
        IdxRadR2dRegionSrcXDetClose = find(radR2dRegionSrcXDet < 40) ;

        clippingRegionSrcXDetByTimeSteps = reshape(clippingRegion,...
            size(clippingRegion,1)*size(clippingRegion,2),size(clippingRegion,3)) ;

        clippingRegionSrcXDetByTimeStepsClose = ...
            zeros(size(clippingRegion,1)*size(clippingRegion,2),...
            size(clippingRegion,3)) ;
        clippingRegionSrcXDetByTimeStepsClose(...
            IdxRadR2dRegionSrcXDetClose,:) = ...
            clippingRegionSrcXDetByTimeSteps(IdxRadR2dRegionSrcXDetClose,:);
        clippingRegionClose = reshape(...
            clippingRegionSrcXDetByTimeStepsClose,size(clippingRegion,1),...
            size(clippingRegion,2),size(clippingRegion,3)) ;

        oFClose=sum(clippingRegionClose(:,:,1:2:end),3);
        eFClose=sum(clippingRegionClose(:,:,2:2:end),3);
        NclipEClose=sum(eFClose(:));
        NclipOClose=sum(oFClose(:));
        

        oF=sum(clippingRegion(:,:,1:2:end),3);
        eF=sum(clippingRegion(:,:,2:2:end),3);
        % figure;subplot(1,2,1);imagesc(oF);title('odd F clipping')
        % xlabel('detectors');ylabel('sources')
        % subplot(1,2,2);imagesc(eF);
        % xlabel('detectors');ylabel('sources');title('even F clipping')

        NclipE=sum(eF(:));
        NclipO=sum(oF(:));

        if NclipEClose~=NclipOClose % cool to replace clipped hot
            if NclipOClose>NclipEClose % 1 is odd is hot, 0 if even is hot
                dHot=MRegion(:,:,:,1:2:end);
                dCool=MRegion(:,:,:,2:2:end);
                % Get subset of clippingRegions since the cool clipped
                % S/D pairs are still clipped
                clippingRegionCool=clippingRegion(:,:,2:2:end);
                fsKeep=fs(2:2:end);
                repM=oF(:)>1;
%                 clippingRegion=eF>0;
                disp(['<<<Interpolating odd/even bright/dim data for Region ',num2str(j+1)])
                info.system.framerate=info.framerate;
                dHrs=resample_tts(dHot,info,2*info.framerate);
                dHot=dHrs(:,:,:,2:2:end);
                dHot(dHot<0)=eps;
                clear dHrs
            else
                dHot=MRegion(:,:,:,2:2:end);
                dCool=MRegion(:,:,:,1:2:end);
                % Get subset of clippingRegions since the cool clipped
                % S/D pairs are still clipped
                clippingRegionCool=clippingRegion(:,:,1:2:end); 
                fsKeep=fs(2:2:end);
                repM=eF(:)>1;
%                 clippingRegion=oF>0;
                disp(['<<<Interpolating even/odd bright/dim data for Region ',num2str(j+1)])
                info.system.framerate=info.framerate;
                dCrs=resample_tts(dCool,info,2*info.framerate);
                dCool=dCrs(:,:,:,2:2:end);
                dCool(dCool<0)=eps;
                clear dCrs
            end
            
            % scale dark data
            for wl=1:Nwl
                dCool(:,:,wl,:)=dCool(:,:,wl,:).*DimScale(wl);
            end
            
            % put dark data in with bright data
            dHot=reshape(dHot,size(dHot,1)*size(dHot,2),2,[]);
            dCool=reshape(dCool,size(dCool,1)*size(dCool,2),2,[]);
            dHot(repM,:,1:Nf_new)=dCool(repM,:,1:Nf_new);


        else % no need for replacement
            M1=reshape(MRegion,Ns*Nd,2,[]);
            keepNN2=find(((Rad.info.pairs.NN==2)==1).*((Rad.info.pairs.WL==1)==1));
            oNN2=mean(mean(mean(M1(keepNN2,:,1:2:end),3),2),1); % odd
            eNN2=mean(mean(mean(M1(keepNN2,:,2:2:end),3),2),1); % even
            if oNN2>eNN2
                dHot=MRegion(:,:,:,1:2:end);
                % Get subset of clippingRegions since the cool clipped
                % S/D pairs are still clipped
                clippingRegionCool=clippingRegion(:,:,2:2:end); 
                fsKeep=fs(1:2:end);
            else
                dHot=MRegion(:,:,:,2:2:end);
                % Get subset of clippingRegions since the cool clipped
                % S/D pairs are still clipped
                clippingRegionCool=clippingRegion(:,:,1:2:end); 
                fsKeep=fs(2:2:end);
            end
        end
        MRegion=reshape(dHot,size(MRegion,1),size(MRegion,2),2,[]);
        if ( j == 0 )
            MBestPassAllRegions = zeros(size(M,1),size(M,2),size(M,3),size(MRegion,4)) ;
            clippingAllRegions = zeros(size(M,1),size(M,2),size(MRegion,4)) ;
        end
        MBestPassAllRegions(RegSrcStartIdx:RegSrcEndIdx,:,:,:) = MRegion ;
        clippingAllRegions(RegSrcStartIdx:RegSrcEndIdx,:,:) = clippingRegionCool ;
    end
    end
    M = MBestPassAllRegions ;
    clipping = clippingAllRegions ;
    info.framerate=info.framerate/2; % % OR UP-SAMPLE????
    info.nframe=size(M,4);
    fs=fsKeep;
end


%% Load in aux channels and stim synch
disp(['<<< Loading Stim Synch data'])
if info.naux
    aux_chan=[(FSchan-1):-1:(FSchan-info.naux)];
    
    % Load aux channels
    for n=1:length(aux_chan)
        % aux file will be named DATE-run###-ch##.raw
        [rawFilename,rawFilenameError] = fnbase2rawFileName(...
            [dateA,'-run',num2str(info.run,'%03g'),'-ch'],...
            aux_chan(n), ext ) ;
        fid=fopen(rawFilename) ;
        
        if fid==-1
            error(['** loaddata: No .raw file was found for aux ',...
                num2str(n),' (channel ',num2str(num2str(aux_chan(n) )),') **'])
        end
        
        aux.(['aux',num2str(n)])=fread(fid,'float32');
        fclose(fid);
    end
    
    % Stim synch
    synch=zeros(2,info.nframe-1);
    ti=fs(1:info.nframe-1);
    tIdx=bsxfun(@plus,ti,[0:floor(info.framesize)-1])';
    s1=reshape(aux.aux1(tIdx(:)),floor(info.framesize),[])';
    synch(1,:)=squeeze(std(s1,[],2));
    fsynch=abs(fft(s1,info.framesize,2));
    [~,synch(2,:)]=max(fsynch(:,1:floor(info.framesize/2)),[],2);    
    
    [info.synchpts,info.synchtype]=findsynch(synch);
    PulseTypes=unique(info.synchtype,'stable');
    for j=1:numel(PulseTypes)
        info.(['Pulse_',num2str(j)])=find(info.synchtype==PulseTypes(j));
    end
    info.synch=synch;
    info.fs=ti;
    info.raw_synch=aux.aux1;
    
else
    synch=[];
end


%% Adjust outputs if not IQ decoding
if ~nnz(I)
    I=[];
    Q=[];
end

%% Get grid and radius files into structure

if exist(['Pad_',info.pad,'.mat'])
    foo=load(['Pad_',info.pad]);
    info.Pad_info=foo.info;
else
    
    foo=load(['grid_',info.pad]);
    if isfield(foo,'grid')
        info.grid=foo.grid;
    else
        info.grid=foo;
    end
    
    foo=load(['radius_',info.pad]);
    if isfield(foo,'Rad')
        info.Rad=foo.Rad;
    else
        info.Rad=foo;
    end
end

if params.ND2 % Output data in newer format 181217)
    % data is source-detector-wavelength-timepoints
    M=reshape(M,[],size(M,4));
    clipping=reshape(clipping,[],size(clipping,3));
    info = converter_info(info, 'ND1 to ND2');
    clipping=cat(1,clipping,clipping);
    clipped=sum(clipping,2);
        info.MEAS= struct;%table(clipped(:)>0,'VariableNames', {'Clipped'});
        info.MEAS.Clipped = clipped(:)>0;
end


%Store all units
if ~isfield(info, 'units')
    info.units = struct;
    info.units.dpos2 = 'mm';
    info.units.dpos3 = 'mm';
    info.units.spos3 = 'mm';
    info.units.spos2 = 'mm';
    info.units.lambda = 'nm';
    info.units.Mod = 'MHz';
    info.units.r2d = 'mm';
    info.units.r3d = 'mm';
    info.units.framerate = 'Hz';
    info.units.init_framerate = 'Hz';
    info.units.mua = 'mm^-1';
    info.units.musp = 'mm^-1';
    info.units.cortex_mu_a = 'mm^-1';
    info.units.cortex_Hb = 'micro_mol';
end

if ( ~params.noMatFile )
save([datafile2save,'.mat'],'M','info','I','Q','clipping','-v7.3')
end
end

