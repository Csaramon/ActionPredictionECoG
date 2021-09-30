%% download file from server
server_ip = '192.87.11.69';
base_path = '/Users/qinchaoyi/Desktop/ActionPrediction/';
server_path = '/data00/Chaoyi/ActionPrediction/';


allsub = {'Patient1','Patient2','Patient3','Patient4','Patient6', ...
    'Patient8','Patient9','Patient11','Patient12','Patient13'}; % all sub

for isub = 1:numel(allsub)
    
    fprintf(['Currently processing: sub' num2str(isub) ' ' allsub{isub} '\n']);
    subpath = [base_path allsub{isub} filesep];
    
    %% download file
    destineypath = [subpath filesep 'Analysis' filesep];
    if ~exist(destineypath,'file')
        mkdir(destineypath)
    end
    downloadfile1 =  [server_path allsub{isub} filesep 'Analysis' filesep allsub{isub} 'PowDiff.mat'];
    unix(['expect expscpDownload.exp ' server_ip ' qin Sblqcy! ' downloadfile1 ' ' destineypath])
    
%     downloadfile2 =  [server_path allsub{isub} filesep 'mri' filesep 'T1.mgz'];
%     unix(['expect expscpDownload.exp 159.226.113.181 cyqin wlabqcy! ' downloadfile2 ' ' destineypath])
%     
%     downloadfile3 =  [server_path allsub{isub} filesep 'brain3D' filesep 'stdspace/elec_rois_individual.nii'];
%     unix(['expect expscpDownload.exp 159.226.113.181 cyqin wlabqcy! ' downloadfile3 ' ' destineypath])
%     
    
    
    
end