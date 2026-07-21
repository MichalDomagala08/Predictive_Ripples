function sourceLoc = findSpikeTimes_reviewed(X_raw,peak_win, z_thresh, cts_thresh,spikeTimeAround,varargin)
%%% Spike Detection algorithm - adapted from Diamond et al., 2010 

%% parse inputs
inp_pars = inputParser;

% specific to this code as a ascale to amplitue
def_amp_scale=3;
addParameter(inp_pars,'amp_scale',def_amp_scale,@(x) isscalar(x));

% For passing in time series. Peaks Characteristics (add default values)
addParameter(inp_pars,'maxNegPeakWidth',50); 
addParameter(inp_pars,'maxPosPeakWidth',Inf); 
addParameter(inp_pars,'peakSeparation',0); 
addParameter(inp_pars,'trackPeaks',false); 
addParameter(inp_pars,'maxPeakHeight',Inf); 

parse(inp_pars,varargin{:}) % parse to get 
amp_scale=inp_pars.Results.amp_scale;

maxNegPeakWidth = inp_pars.Results.maxNegPeakWidth;
maxPosPeakWidth = inp_pars.Results.maxPosPeakWidth;
trackPeaks = inp_pars.Results.trackPeaks; 
peakSeparation = inp_pars.Results.peakSeparation; 
maxPeakHeight = inp_pars.Results.maxPeakHeight; 

if trackPeaks %Starting from Positive PEak
    posPeakSeparation = peakSeparation; negPeakSeparation = 0; 
    posPeakHeight = z_thresh; negPeakHeight = 0; 
else %Starting from negative trough
    negPeakSeparation = peakSeparation; posPeakSeparation = 0; 
    posPeakHeight = 0; negPeakHeight = z_thresh; 
end

%% load data and set up parameters/arrays

if ~isstruct(X_raw); error('Check input arguments. X_raw should be a struct now, with Fs field.'); end
Fs = X_raw.Fs;
chanNames = X_raw.chanNames;
X_raw = X_raw.rawTs;

dims=size(X_raw);

hp=floor(peak_win*Fs); % window around peak 0.1 in seconds * 500Hz - 50

%% find peaks based on polarity
warning('off','signal:findpeaks:largeMinPeakHeight');

fullRaster = sparse(dims(1),dims(2));
waveformsMaster = cell(1,dims(2)); 

%for every channel
for kk=1:length(chanNames)
    
    currentRaster = sparse(dims(1),1);
    waveforms = cell(size(currentRaster)); 

    %% Peaks and troughs 
    
    %Find firstly Positive Peak (If it is established to do so, that
    %exceeded posPeakHeight threshold and have duration of utmost
    %posPeakSeparation duration

    % Find Peaks Description:
    %   + Prominence:Realtive importance of a peak:
            %  Marker on a peak, extending horizontal ine  to the point when it corsses the signal in a gighest peak,
            %  or reaches left or right end of the signal 
            %Then we fund the minimum signal in both those intervals  - So
            %it is wysokość względna 
%        Here, Promiennce is at Least 4 std !!!!! bREMEMBER that signal is
%        ZSCORED beforehand
%       + 

    xCurrent = X_raw(:,kk);
    [pHeight,pInd]=findpeaks(xCurrent,'MinPeakprominence',z_thresh,'MinPeakHeight',posPeakHeight,'MaxPeakWidth',maxPosPeakWidth,'minPeakDistance',posPeakSeparation);
    %%% Prominence: 4 std; min Positive Peak Height - 0 ( All) , maxwidh -Indefinite; Min Peak Distance - 0 (All) 
    % So it detects ALL sufficiently High Postiive Peaks) 

    % Checkingig if the threshold is exceeded By maximum (Infinte
    % (unattainable)) If so, then zero them! --> To mitigate Bad Signal
    isBad = pHeight > maxPeakHeight; 
    pHeight(isBad) = []; 
    pInd(isBad) = []; 
    

    %Find firstly Negative Peak (If it is established to do so), that
    %exceeded negPeakHeight threshold and have duration of utmost
    %negPeakSeparation duration
    [nHeight,nInd]=findpeaks(-xCurrent,'MinPeakProminence',z_thresh,'MinPeakHeight',negPeakHeight,'MaxPeakWidth',maxNegPeakWidth,'minPeakDistance',negPeakSeparation);
    % reverse xCurrent, Prominence 4 std, Min height - 4 std, maxWidth = 50
    % (???- Default - 50 std is a lot); No separation
    isBad = false(size(nHeight));
    nHeight(isBad) = [];
    nInd(isBad) = [];


    %% Matching
        
    maxNumel = 1000000; 
    bufferSize = min(100, ceil(maxNumel / length(pInd))); %How many Peaks are computed together?
    
    nIndBuffer = buffer(nInd,bufferSize);    %Buffering signal Indicies into smaller szize (Negative)
    Nh2Buffer = buffer(nHeight,bufferSize);  % Buffering Signal Height
    
    Nh2Buffer(nIndBuffer == 0) = NaN; %Remove all zeros for NaNs
    nIndBuffer(nIndBuffer == 0) = NaN;
    
    %%% Trying to find matching Positive/Negative Peak in at most 100ms
    %%% window from our current spike 
    for ii = 1:size(nIndBuffer,2)
        
        matchingMatrix = pInd - nIndBuffer(:,ii)'; % make Matrix Matching All Negative and Positive Peaks Together
        
        % inds = matchingMatrix >= -hp & matchingMatrix < 0;
        % The above -- for up-deflection to come before down-deflection.
        inds = matchingMatrix >= -hp & matchingMatrix <= hp;
        
        [pFind,nFind] = find(inds); % Find height and width coordinates of matched Peaks
        
        heightMatrix = false(size(pFind));
        for jj = 1:length(pFind) %for every Height of Positive and NEgative togerher - check whether they exceed amp_scale (1) * threshold 4  
            heightMatrix(jj) = pHeight(pFind(jj)) + Nh2Buffer(nFind(jj),ii) >= amp_scale * z_thresh;
        end
        
        nFind = nFind(heightMatrix); 
        pFind = pFind(heightMatrix); 
                
        % Waveform
        if trackPeaks %Start from Positive
            for jj = 1:length(pFind)
                ts = nan(1,2 * hp + 1);
                % tsBb = ts; 
                
                startInd = max(pInd(pFind(jj)) - hp,1);
                startBuffer = max(-(pInd(pFind(jj)) - hp) + 2,1);
                
                endInd = min(pInd(pFind(jj)) + hp,dims(1));
                endBuffer = min(hp-(pInd(pFind(jj))-dims(1)) + 1,2 * hp + 1);
                % ts = xCurrent(pInd(pFind(jj))- hp:pInd(pFind(jj)) + hp);
                ts(startBuffer:endBuffer) = xCurrent(startInd:endInd);
                waveforms{pInd(pFind(jj))} = ts;
            end
        else %Sstart from Negative 
            for jj = 1:length(nFind)
                
                ts = nan(1,2 * hp + 1);
                % I'll start by filling this with NaNs, for the rare event that
                % we have spikes at the very beginning of the time series.
                % tsBb = ts; 
                
                
                startInd = max(nIndBuffer(nFind(jj),ii) - hp,1);
                startBuffer = max(-(nIndBuffer(nFind(jj),ii) - hp) + 2,1);
                
                endInd = min(nIndBuffer(nFind(jj),ii) + hp,dims(1));
                endBuffer = min(hp-(nIndBuffer(nFind(jj),ii)-dims(1)) + 1,2 * hp + 1);
                % ts = xCurrent(nIndBuffer(nFind(jj),ii)- hp:nIndBuffer(nFind(jj),ii) + hp);
                ts(startBuffer:endBuffer) = xCurrent(startInd:endInd);
                
                %             clf; plot(nIndBuffer(nFind(jj),ii)- hp:nIndBuffer(nFind(jj),ii) + hp,ts); hold on
                %             plot(nIndBuffer(nFind(jj),ii),xCurrent(nIndBuffer(nFind(jj),ii)),'ro')
                %             plot(pInd(pFind(jj)),xCurrent(pInd(pFind(jj))),'bo');
                % pause
                
                waveforms{nIndBuffer(nFind(jj),ii)} = ts;
                                
            end
        end
        
        if trackPeaks; currentRaster(pInd(pFind)) = true; % To retain peaks 
        else; currentRaster(nIndBuffer(nFind,ii)) = true;
        end
        
    end
    waveforms(~currentRaster) = []; 
    waveforms = cell2mat(waveforms); 
    waveformsMaster{kk} = waveforms; 
    
    fullRaster(:,kk) = currentRaster;
    
    
end



%% Narrow by counts 

badCounts = sum(fullRaster) < cts_thresh; 
fullRaster(:,badCounts) = false; 

%% Spike Time Computation %%% CHECK!!!
for chan = 1:size(fullRaster,2)
    spikeTime{chan} = [find(fullRaster(:,chan)==1)-spikeTimeAround,find(fullRaster(:,chan)==1)+spikeTimeAround];
end
%% Pack up 

sourceLoc.rasters = fullRaster;
sourceLoc.waveforms = waveformsMaster;
sourceLoc.spikeTime = spikeTime;

end

