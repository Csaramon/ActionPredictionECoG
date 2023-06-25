function shufffledData = shuffleData(inputData,Fs,numSec)

% inputData should be trials*timepoints
if nargin <3
    numSec = ceil((size(inputData,2))/Fs);
end

shufffledData = [];
for itrl = 1:size(inputData,1)
    
    shfTrlData = [];
    trlData = inputData(itrl,:)';
    % Choose numSec random 'cut' positions
    dpsplit = ceil(size(trlData,1).*rand(numSec,1));
    % Arrange these in ascending order
    dpsplit = sort (dpsplit);
    
    start(1)=1;
    start(2:numSec)=dpsplit(1:numSec-1);
    ending(1:numSec-1)=dpsplit(1:numSec-1)-1;
    ending(numSec) =  size(trlData,1);
    
    order = randperm(numSec);
    
    for c = 1:numSec
        
        %shuffle the signal
        shfTrlData = [shfTrlData; trlData(start(order(c)):ending(order(c)),:)];
        
    end
    
    
    shufffledData(itrl,:) = shfTrlData'; 
    
    
end




end % function end