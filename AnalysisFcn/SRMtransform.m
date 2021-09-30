function [srm] = SRMtransform(origdata,iteration,feature)

% input:
% origdata should be 3 dimensional data (electrodes*timepoints*subjects)
% use Nan to represent missing electrodes 

% output:
% srm is also 3 dimensional data (subjects*feastures*timepoints)

if nargin < 3
    feature = min(sum(~isnan(squeeze(origdata(:,1,:))),1));
    if nargin < 2
        iteration = 10;
    end
end

if feature > min(sum(~isnan(squeeze(origdata(:,1,:))),1))
    warning(['Number of features should not exceed the minimum valid number of electrodes!' ...
    'So it was automatically adjusted.'])
    feature = min(sum(~isnan(squeeze(origdata(:,1,:))),1));
end
    
save('SRMdatatmp','origdata','-v6')

[s,w] = unix(['python3 SRMtransform.py SRMdatatmp.mat ' num2str(iteration) ' ' num2str(feature)]);
if s
    error(w);
end
delete('SRMdatatmp.mat')

tmpFile = load('srmData.mat');
srm = tmpFile.('srmdata');


