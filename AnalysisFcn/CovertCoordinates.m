[CTfilename, CTpathname, ~] = uigetfile(['/Users/qinchaoyi/Desktop/*.nii'],'Choose CT scan to register');
[MRIfilename, MRIpathname, ~] = uigetfile(['/Users/qinchaoyi/Desktop/*.nii'],'Choose template MRI scan');

VCT = spm_vol([CTpathname CTfilename]);
VMRI = spm_vol([MRIpathname MRIfilename]);

[ROIpathname] = uigetdir('/Users/qinchaoyi/Desktop', 'Pick the directory containing the ROI labels');

ROIFiles = dir([ROIpathname filesep '*.label']);

spm_figure('Create','Graphics','coregistration',1);
flags.graphics = 1;
x = spm_coreg(VMRI,VCT,flags);
Tmat = VCT.mat\spm_matrix(x(:)')*VMRI.mat;


for iroi = 1:numel(ROIFiles)
CTCoor = read_label([],[ROIFiles(iroi).folder filesep ROIFiles(iroi).name]);
   CTCoor = CTCoor(:,2:4);

MRICoor = VMRI.mat*(Tmat\(VCT.mat\[CTCoor ones(size(CTCoor,1),1)]'));
MRICoor = MRICoor(1:3,:)';


end
