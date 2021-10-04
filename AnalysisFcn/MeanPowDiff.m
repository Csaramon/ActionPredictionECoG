% Calculate grand average power difference with permutation test
%
 
% input data time series and sampling rate
% allTS = rerefData.trial{1};
% fs = rerefData.fsample;
 
% number of permutations
Nshf = 2000;
 
%% reconcantenate time series for the two views in two conditions
 
% normalize across frequency points
mu =  mean(allTS,3,'omitnan');
allTS = allTS./repmat(mu,1,1,size(allTS,3));
allTS = squeeze(mean(allTS,2,'omitnan'));
 
% predefine variable for the reconcantenate time series
tsI = [];
tsS = [];
 
%%%%-------- reconcatenate data movie-wise and condition-wise --------%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % extract each movie's start time point (with first camera change excluded)
% nc = 1;
% while nc < size(camInfo,1)
%     
%     if camInfo{nc,3}==camInfo{nc+1,3}
%         camInfo(nc+1,:) = [];
%     else
%         nc = nc+1;
%     end
%     
% end
% 
% % all the unique movies
% uniMovie = unique(cell2mat(camInfo(:,3)));
% 
% % reconcatenate data movie-wise and condition-wise
% for im = uniMovie'
%     
%     fInd = find(cell2mat(camInfo(:,3))==im);
%     for ii = fInd'
%         % extract data of each movie
%         tmpTS = allTS(:,camInfo{ii,4}:camInfo{ii,4}+fs*(camInfo{ii,7}-camInfo{ii,5}));
%         
%         
%         % concatenate data according to condition
%         if camInfo{ii,1}==0
%             tsI=cat(2,tsI,mean(tmpTS,2,'omitnan'));
%         elseif camInfo{ii,1}==1
%             tsS=cat(2,tsS,mean(tmpTS,2,'omitnan'));
%         end
%         
%         
%     end
%     
% end
 
%%%%-------- reconcatenate data movie clip-wise and condition-wise --------%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for im = 1:size(camInfo,1)
    
    
    % extract data of each movie
    tmpTS = allTS(:,camInfo{im,4}:camInfo{im,4}+round(fs*camInfo{im,6}));
    
    
    % concatenate data according to condition
    if camInfo{im,1}==0
        tsI=cat(2,tsI,mean(tmpTS,2,'omitnan'));
    elseif camInfo{im,1}==1
        tsS=cat(2,tsS,mean(tmpTS,2,'omitnan'));
    end
    
    
    
end
 
%% calculate mean power difference with permutation test
if strcmp(calculate,'MPD')
    % calculate the difference
    
    [~,p,~,tt]=ttest2(tsI',tsS');
    PD_I_S = tt.tstat';
    
    % predefine variable for the null distribution of correlation matrix
    
    PD_shf = zeros(size(tsI,1),Nshf);
    
    strlen = 0;
    % calculate the null distribution of power difference
    for ishf = 1:Nshf
        
        s = ['Permutation progress:' num2str(ishf) '/' num2str(Nshf)];
        strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
        strlen = strlentmp - strlen;
        
        randInd = randperm(size(tsI,2)+size(tsS,2));
        tsTmp = cat(2,tsI,tsS);
        tsIrand = tsTmp(:,randInd(1:size(tsI,2)));
        tsSrand = tsTmp(:,randInd(size(tsI,2)+1:end));
        [~,~,~,ttrand]=ttest2(tsIrand',tsSrand');
        PD_shf(:,ishf) =   ttrand.tstat;
        
        
    end
    
    % calculate P value
    
    PD_shfsort = sort(PD_shf,2);
    [~,im] = min(abs(repmat(PD_I_S,1,Nshf)-PD_shfsort),[],2);
    PD_I_S(:,2) = 1-im/Nshf;
end
