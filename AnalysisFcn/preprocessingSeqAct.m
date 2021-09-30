%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%---- Sequential Action Movie Task ----%%%%%%%%
%%%%%%%%------------ Preprocessing ---------%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Directory Configurations

clear;clc;
% specify path for fieldtrip toolbox
addpath('/Users/qinchaoyi/Documents/MATLAB/toolbox/fieldtrip-20210418/')
ft_defaults

% specify datapath
dataPath = '/Users/qinchaoyi/Desktop/ActionPrediction/';

rerefMethod = 'LAR';

patientID = 'Patient2';  % 1,2,3,4,6,8,9
mkdir([dataPath patientID],'Analysis');

% specify paths for recording data and stimuli info
rawdataDir = [dataPath patientID filesep 'RawData'];
stimDir = [dataPath patientID filesep 'Stimuli'];
elecDir = [dataPath patientID filesep 'FsRecon/brain3D'];
resultDir = [dataPath patientID filesep 'Analysis'];

% initialize data matrix
rawdataFiles = dir([rawdataDir filesep '*.mat']);
rawdataComb = [];
triggerComb = [];

switch patientID
    
    case 'Patient1'
        % manually correct the channel label and trigger channel
        for fileNum = 1:size(rawdataFiles,1)
            
            load([rawdataFiles(fileNum).folder filesep rawdataFiles(fileNum).name],'hdr','raw','trigger_photodiode')
            rawdataComb(:,end+1:end+size(raw,2)) = raw;
            triggerComb(:,end+1:end+size(trigger_photodiode,2)) = trigger_photodiode;
            clear raw trigger_photodiode
        end
        
        % load electrode coordinates
        load([elecDir filesep 'ft_elec.mat'])
        % remove redundant channels and correct labels
        chan2remove = [59]; % reference channel
        rawdataComb(chan2remove,:) = [];
        hdr.nSamples = size(rawdataComb,2);
        hdr.orig.epochdef(1,2) = size(rawdataComb,2);
        hdr.nChans = size(rawdataComb,1);
        hdr.label = elec.label;
        hdr.chantype = elec.chantype;
        hdr.chanunit = elec.chanunit;
        
    case 'Patient2'
        % manually correct the channel label and trigger channel
        for fileNum = 1:size(rawdataFiles,1)
            
            load([rawdataFiles(fileNum).folder filesep rawdataFiles(fileNum).name],'hdr','raw')
            trigger = raw(63,:);
            rawdataComb(:,end+1:end+size(raw,2)) = raw;
            triggerComb(:,end+1:end+size(trigger,2)) = trigger;
            clear raw trigger
        end
        
        % load electrode coordinates
        load([elecDir filesep 'ft_elec.mat'])
        % remove redundant channels and correct labels
        chan2remove = [63,126];
        rawdataComb(chan2remove,:) = [];
        hdr.nSamples = size(rawdataComb,2);
        hdr.orig.epochdef(1,2) = size(rawdataComb,2);
        hdr.nChans = size(rawdataComb,1);
        hdr.label = elec.label;
        hdr.chantype = elec.chantype;
        hdr.chanunit = elec.chanunit;
        
    case 'Patient3'
        % manually correct the channel label and trigger channel
        for fileNum = 1:size(rawdataFiles,1)
            
            load([rawdataFiles(fileNum).folder filesep rawdataFiles(fileNum).name],'hdr','raw','trigger_photodiode')
            rawdataComb(:,end+1:end+size(raw,2)) = raw;
            triggerComb(:,end+1:end+size(trigger_photodiode,2)) = trigger_photodiode;
            clear raw trigger_photodiode
        end
        
        % load electrode coordinates
        load([elecDir filesep 'ft_elec.mat'])
        % remove redundant channels and correct labels
        chan2remove = [];
        rawdataComb(chan2remove,:) = [];
        hdr.nSamples = size(rawdataComb,2);
        hdr.orig.epochdef(1,2) = size(rawdataComb,2);
        hdr.nChans = size(rawdataComb,1);
        hdr.label = elec.label;
        hdr.chantype = elec.chantype;
        hdr.chanunit = elec.chanunit;
        
    case 'Patient4'
        % manually correct the channel label and trigger channel
        for fileNum = 1:size(rawdataFiles,1)
            
            load([rawdataFiles(fileNum).folder filesep rawdataFiles(fileNum).name],'hdr','raw','trigger_photodiode')
            rawdataComb(:,end+1:end+size(raw,2)) = raw;
            triggerComb(:,end+1:end+size(trigger_photodiode,2)) = trigger_photodiode;
            clear raw trigger_photodiode
        end
        
        % load electrode coordinates
        load([elecDir filesep 'ft_elec.mat'])
        % remove redundant channels and correct labels
        chan2remove = [125];
        rawdataComb(chan2remove,:) = [];
        hdr.nSamples = size(rawdataComb,2);
        hdr.orig.epochdef(1,2) = size(rawdataComb,2);
        hdr.nChans = size(rawdataComb,1);
        hdr.label = elec.label;
        hdr.chantype = elec.chantype;
        hdr.chanunit = elec.chanunit;
        
    case 'Patient6'
        % manually correct the channel label and trigger channel
        for fileNum = 1:size(rawdataFiles,1)
            
            load([rawdataFiles(fileNum).folder filesep rawdataFiles(fileNum).name],'hdr','raw','trigger')
            
            rawdataComb(:,end+1:end+size(raw,2)) = raw;
            triggerComb(:,end+1:end+size(trigger,2)) = trigger;
            clear raw trigger
        end
        
        % load electrode coordinates
        load([elecDir filesep 'ft_elec.mat'])
        % remove redundant channels and correct labels
        chan2remove = [121,122];
        rawdataComb(chan2remove,:) = [];
        hdr.nSamples = size(rawdataComb,2);
        hdr.orig.epochdef(1,2) = size(rawdataComb,2);
        hdr.nChans = size(rawdataComb,1);
        hdr.label = elec.label;
        hdr.chantype = elec.chantype;
        hdr.chanunit = elec.chanunit;
        
    case 'Patient8'
        % manually correct the channel label and trigger channel
        for fileNum = 1:size(rawdataFiles,1)
            
            load([rawdataFiles(fileNum).folder filesep rawdataFiles(fileNum).name],'hdr','raw','trigger')
            
            rawdataComb(:,end+1:end+size(raw,2)) = raw;
            triggerComb(:,end+1:end+size(trigger,2)) = trigger;
            clear raw trigger
        end
        
        % load electrode coordinates
        load([elecDir filesep 'ft_elec.mat'])
        % remove redundant channels and correct labels
        chan2remove = 107;
        rawdataComb(chan2remove,:) = [];
        hdr.nSamples = size(rawdataComb,2);
        hdr.orig.epochdef(1,2) = size(rawdataComb,2);
        hdr.nChans = size(rawdataComb,1);
        hdr.label = elec.label;
        hdr.chantype = elec.chantype;
        hdr.chanunit = elec.chanunit;
        
    case 'Patient9'
        % manually correct the channel label and trigger channel
        for fileNum = 1:size(rawdataFiles,1)
            
            load([rawdataFiles(fileNum).folder filesep rawdataFiles(fileNum).name],'hdr','raw','trigger')
            
            rawdataComb(:,end+1:end+size(raw,2)) = raw;
            triggerComb(:,end+1:end+size(trigger,2)) = trigger;
            clear raw trigger
        end
        
        % load electrode coordinates
        load([elecDir filesep 'ft_elec.mat'])
        % remove redundant channels and correct labels
        chan2remove = [45 46 51 52 57 58 63];
        rawdataComb(chan2remove,:) = [];
        hdr.nSamples = size(rawdataComb,2);
        hdr.orig.epochdef(1,2) = size(rawdataComb,2);
        hdr.nChans = size(rawdataComb,1);
        hdr.label = elec.label;
        hdr.chantype = elec.chantype;
        hdr.chanunit = elec.chanunit;
        
end

raw = rawdataComb;
trigger = triggerComb;
clear rawdataComb triggerComb

%% Load Movie Information

movieSeq = cell(8,size(rawdataFiles,1));
for fileNum = 1:size(rawdataFiles,1)
    switch rawdataFiles(fileNum).name(end-5)
        case 'A'
            tmpFile = fopen([stimDir filesep 'StimulusList_A.txt'],'r');
            tmpMovieList = textscan(tmpFile,'%s');
            movieSeq(:,fileNum) = tmpMovieList{1,1};
            fclose(tmpFile);
            clear tmpFile tmpMovieList
        case 'B'
            tmpFile = fopen([stimDir filesep 'StimulusList_B.txt'],'r');
            tmpMovieList = textscan(tmpFile,'%s');
            movieSeq(:,fileNum) = tmpMovieList{1,1};
            fclose(tmpFile);
            clear tmpFile tmpMovieList
        case 'C'
            tmpFile = fopen([stimDir filesep 'StimulusList_C.txt'],'r');
            tmpMovieList = textscan(tmpFile,'%s');
            movieSeq(:,fileNum) = tmpMovieList{1,1};
            fclose(tmpFile);
            clear tmpFile tmpMovieList
    end
end

movieSeq = movieSeq(:);
numMovies = length(movieSeq);

uniMovieList = unique(movieSeq);
for movie = 1:length(uniMovieList)
    uniMovieList{movie,2} = movie;
end


%% Setup Trigger Threshold and Detect Trigger Samples

figure
plot(trigger)
[~,trigThresh] = ginput(1);
close

plot(trigger*2)
hold on
tmpTrigger = trigger;
tmpTrigger(tmpTrigger<trigThresh) = 0;
findpeaks(tmpTrigger);

% detect trigger samples
[~,movieSamples] = findpeaks(tmpTrigger,'MinPeakDistance',100);

if length(movieSamples)==numMovies
    disp('The number of triggers matches with the number of movies.');
else
    warning(' The number of triggers does not match with the movies.');
end

clear triggerComb tmpTrigger trigThresh
%% Extract Camera Change Information

cams2Rej = 1;

camIndex = 1;
for movie = 1:numMovies
    tmpStimInfo = readtable([stimDir filesep movieSeq{movie,1} '.csv']);
    
    for cam = cams2Rej+1:size(tmpStimInfo,1)
        if  strcmp(movieSeq{movie,1}(1,1),'M')
            camInfo{camIndex,1} = 0;
        elseif strcmp(movieSeq{movie,1}(1,1),'S')
            camInfo{camIndex,1} = 1;
        end
        camInfo{camIndex,2} = movieSeq{movie,1};
        camInfo{camIndex,3} = uniMovieList{strcmp(movieSeq{movie,1},uniMovieList(:,1)),2};
        camInfo{camIndex,4} = movieSamples(1,movie)+(tmpStimInfo{cam,12}.*hdr.Fs);
        camInfo{camIndex,5} = tmpStimInfo{cam,12}-tmpStimInfo{cam-1,12};
        if cam+1>size(tmpStimInfo,1)
            camInfo{camIndex,6} = tmpStimInfo{cam,4}-tmpStimInfo{cam,12};
        else
            camInfo{camIndex,6} = tmpStimInfo{cam+1,12}-tmpStimInfo{cam,12};
        end
        camInfo{camIndex,7} = tmpStimInfo{cam,4};
        
        camIndex = camIndex+1;
    end
    clear tmpStimInfo
end

clear cams2Rej camIndex movie cam


%% Transform Data into FieldTrip Format

data = [];
data.hdr = hdr;
data.trial{1} = raw;
data.fsample = hdr.Fs;
data.label = hdr.label;
data.time{1} = (1:hdr.nSamples)./hdr.Fs; % in seconds
data.sampleinfo = [1 hdr.nSamples];
data.elec = elec;
data.cfg = [];

clear hdr raw trigger

cfg.channel = 'all';
ftData = ft_selectdata(cfg,data);

% save rawdata with trial info
uisave('ftData','camInfo', [resultDir filesep patientID '_rawData.mat']);
% save eventdata
uisave('camInfo', [resultDir filesep patientID '_eventdata.mat']);


%% Seperate Trials

preTrigger = .5;
postTrigger = 1;

numCams = size(camInfo,1);
trl = nan(numCams,8);
for cam = 1:numCams
    trl(cam,1:8) = [camInfo{cam,4}-(preTrigger*ftData.fsample)...
        (camInfo{cam,4}+(postTrigger*ftData.fsample)-1)...
        -preTrigger*ftData.fsample...
        camInfo{cam,1}...
        camInfo{cam,3}...
        camInfo{cam,4}...
        camInfo{cam,5}...
        camInfo{cam,6}];
end

% construct fieldtrip trial data
cfg = [];
cfg.trl = trl;
trlData = ft_redefinetrial(cfg,ftData);

uisave('trlData', [resultDir filesep patientID '_trlData.mat']);

%% inspect spectrum

eeglab;close % use the function of eeglab toolbox
figure;
tmpDat = cell2mat(trlData.trial);
conDat = reshape(tmpDat,numel(trlData.label), numel(trlData.time{1}), []);
spectopo(conDat,size(conDat,2),round(trlData.fsample),'freqrange',[0 250],'percent',50);

%% Remove line noise according to spectum
% preprocess data
cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = 'all';
cfg.lpfilter = 'yes';
cfg.lpfreq = 200;
cfg.padding = 2;
cfg.padtype = 'data';
% remove line noise (band stop filtering)
cfg.bsfilter = 'yes';
cfg.bsfiltord = 3;

% manually specify line noise frequencies and bad channels
switch patientID
    case 'Patient1'
        cfg.bsfreq = [98 102;148 152]; % line noise frequency [48 52;98 102;148 152]
        badChan = [105];
    case 'Patient2'
        cfg.bsfreq = [98 102;148 152]; % line noise frequency [48 52;98 102;148 152]
        badChan = [];
    case {'Patient3','Patient4','Patient8','Patient9'}
        cfg.bsfreq = [48 52;98 102;148 152]; % line noise frequency [48 52;98 102;148 152]
        badChan = [];
    case 'Patient6'
        cfg.bsfreq = [48 52;98 102;148 152]; % line noise frequency [48 52;98 102;148 152]
        badChan = [9 11 109 110];

end
    
% remove line noise (spectrum interpolation)
% cfg.dftfilter = 'yes';
% cfg.dftreplace = 'neighbour';
% cfg.dftfreq = [50 100 150];
% cfg.dftbandwidth = [2 2 2];
% cfg.dftneighbourwidth = [2 3 4];

for i=1:length(trlData.trial)
    trlData.trial{i}(badChan,:) = 0;
end
      
trlData = ft_preprocessing(cfg,trlData);

figure;
tmpDat = cell2mat(trlData.trial);
conDat = reshape(tmpDat,numel(trlData.label), numel(trlData.time{1}), []);
spectopo(conDat,size(conDat,2),round(trlData.fsample),'freqrange',[0 250],'percent',50);

%% Rereference


if strcmpi(rerefMethod,'LAR')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%-------- Local Average Reference (LAR) --------%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rerefRadius = 12; % Local range 12mm
    distances = squareform( pdist(ftData.elec.chanpos) ); %get all pairwise distances between electrodes, organized in a 2D matrix
    
    rerefMap = (distances < rerefRadius)-eye(length(distances)); %get all electrode within radius distance of electrode to reref
    rerefMap(:,badChan) = 0; % will not use bad channels for rereference
    rerefMap(badChan,:) = 0; 
    neighbours = repmat( - sum(rerefMap, 2) , [1,length(rerefMap)]); %count the nuber of neighbouring electrode to electrode to reref
    neighbours(neighbours==0) = -1;
    rerefMap = eye(length(rerefMap)) + rerefMap ./ neighbours; %generate weights coef matrix map for local rereferencing

    cfg = [];
    cfg.channel = ftData.label;
    cfg.montage.labelold = cfg.channel;
    cfg.montage.labelnew = cfg.channel;
    
    cfg.montage.tra = rerefMap;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%-------- Local Average Reference (LAR) --------%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%-------- Common Average Reference (CAR) --------%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch patientID
        case 'Patient6'
            rerefMethod = 'CommonAverageAndBipolar';
            % specify depth electrodes' index
            depthIndex = [101:104;105:108];
            nonDepthIndex = setdiff(1:numel(ftData.label),depthIndex(:));
            traShade = -1/numel(nonDepthIndex)*ones(numel(cfg.channel));
            traShade(depthIndex(:),:) = 0;
            traShade(:,depthIndex(:)) = 0;
            for idepth = 1:size(depthIndex,1)
                temIndex = depthIndex(idepth,1:end-1);
                traShade(temIndex,temIndex+1) = ...
                    traShade(temIndex,temIndex+1) -1.*eye(numel(temIndex));
            end
            
            cfg = [];
            cfg.channel = ftData.label;
            cfg.montage.labelold = cfg.channel;
            cfg.montage.labelnew = cfg.channel;
            
        case 'Patient9'
            rerefMethod = 'CommonAverageAndBipolar';
            % specify depth electrodes' index
            depthIndex = [41:44;45:48;49:52;53:56];
            nonDepthIndex = setdiff(1:numel(ftData.label),depthIndex(:));
            traShade = -1/numel(nonDepthIndex)*ones(numel(cfg.channel));
            traShade(depthIndex(:),:) = 0;
            traShade(:,depthIndex(:)) = 0;
            for idepth = 1:size(depthIndex,1)
                temIndex = depthIndex(idepth,1:end-1);
                traShade(temIndex,temIndex+1) = ...
                    traShade(temIndex,temIndex+1) -1.*eye(numel(temIndex));
            end
            
            cfg = [];
            cfg.channel = ftData.label;
            cfg.montage.labelold = cfg.channel;
            cfg.montage.labelnew = cfg.channel;
            
        otherwise
            rerefMethod = 'CommonAverage';
            cfg = [];
            cfg.channel = ftData.label;
            cfg.montage.labelold = cfg.channel;
            cfg.montage.labelnew = cfg.channel;
            traShade = -1/numel(cfg.channel)*ones(numel(cfg.channel));
            
    end
    
    cfg.montage.tra = eye(numel(cfg.channel))+traShade;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%-------- Common Average Reference(CAR) --------%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

cfg.updatesens = 'yes';
rerefData = ft_preprocessing(cfg, trlData);

uisave('rerefData',[resultDir filesep patientID rerefMethod '_rerefData.mat']);

%% Data Inspection

cfg = [];
cfg.viewmode = 'vertical';
ft_databrowser(cfg, rerefData);

%% Remove the Detected Trials

cfg = [];
cfg.method = 'summary';
rerefData = ft_rejectvisual(cfg,rerefData);


%% Downsample Data
cfg = [];
cfg.resamplefs = 500;
cfg.detrend         = 'no';
[preprocData] = ft_resampledata(cfg, rerefData);

uisave('preprocData', [resultDir filesep patientID rerefMethod '_preprocData.mat']);

