% Calculate inter-view correlation with permutation test
%
% (Each movie was viewed 2 times)


% downsampling 
cfg = [];
cfg.resamplefs = 500;
cfg.detrend         = 'yes';
cfg.demean         = 'yes';
[rerefData] = ft_resampledata(cfg, rerefData);

% extract data time series and sampling rate
allTS = rerefData.trial{1};
fs = rerefData.fsample;

% number of permutations
Nshf = 1000;

%% reconcantenate time series for the two views in two conditions

% predefine variable for the reconcantenate time series
ts1I = [];
ts2I = [];
ts1S = [];
ts2S = [];

% extract each movie's start time point (with first camera change excluded)
nc = 1;
while nc < size(camInfo,1)
    
    if camInfo{nc,3}==camInfo{nc+1,3}
        camInfo(nc+1,:) = [];
    else
        nc = nc+1;
    end
    
end

% all the unique movies
uniMovie = unique(cell2mat(camInfo(:,3)));

% reconcantenate data movie-wise and condition-wise
for im = uniMovie'
    
    fInd1 = find(cell2mat(camInfo(:,3))==im,1,'first');
    if camInfo{fInd1,1}==0
        ts1I=[ts1I,allTS(:,camInfo{fInd1,4}:camInfo{fInd1,4}+fs*(camInfo{fInd1,7}-camInfo{fInd1,5}))];
    elseif camInfo{fInd1,1}==1
        ts1S=[ts1S,allTS(:,camInfo{fInd1,4}:camInfo{fInd1,4}+fs*(camInfo{fInd1,7}-camInfo{fInd1,5}))];
    end
    
    fInd2 = find(cell2mat(camInfo(:,3))==im,1,'last');
    if camInfo{fInd2,1}==0
        ts2I=[ts2I,allTS(:,camInfo{fInd2,4}:camInfo{fInd2,4}+fs*(camInfo{fInd2,7}-camInfo{fInd2,5}))];
    elseif camInfo{fInd1,1}==1
        ts2S=[ts2S,allTS(:,camInfo{fInd2,4}:camInfo{fInd2,4}+fs*(camInfo{fInd2,7}-camInfo{fInd2,5}))];
    end
    
    
end

%% calculate interview correaltion with permutation test

% calculate the correlation matrix
IVC_I = diag(corr(ts1I',ts2I'));
IVC_S = diag(corr(ts1S',ts2S'));

% predefine variable for the null distribution of correlation matrix
IVC_Ishf = zeros(size(ts1I,1),Nshf);
IVC_Sshf = zeros(size(ts1I,1),Nshf);

strlen = 0;
% calculate the null distribution of correlation matrix
for ishf = 1:Nshf
    
        s = ['Permutation progress:' num2str(ishf) '/' num2str(Nshf)];
    strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
    strlen = strlentmp - strlen;
    
    randInd = randi(length(ts1I),1);
    IVC_Ishf(:,ishf) = diag(corr(ts1I(:,[randInd:end,1:randInd-1])',ts2I'));
    
    randInd = randi(length(ts1S),1);
    IVC_Sshf(:,ishf) = diag(corr(ts1S(:,[randInd:end,1:randInd-1])',ts2S'));
       
end

% calculate P value
IVC_Ishfsort = sort(IVC_Ishf,2);
[~,im] = min(abs(repmat(IVC_I,1,Nshf)-IVC_Ishfsort),[],2);
IVC_I(:,2) = 1-im/Nshf;

IVC_Sshfsort = sort(IVC_Sshf,2);
[~,im] = min(abs(repmat(IVC_I,1,Nshf)-IVC_Sshfsort),[],2);
IVC_S(:,2) = 1-im/Nshf;