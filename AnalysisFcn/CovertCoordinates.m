function CovertCoordinates
[CTfilename, CTpathname, ~] = uigetfile(['/Users/qinchaoyi/Desktop/*.nii'],'Choose CT scan to register');
[MRIfilename, MRIpathname, ~] = uigetfile(['/Users/qinchaoyi/Desktop/*.nii'],'Choose template MRI scan');

VCT = spm_vol([CTpathname CTfilename]);
CTVOL = spm_read_vols(VCT);
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
    
    calCTCoor = centroidCalibration(CTCoor,CTVOL,VCT.mat);
    
    
    MRICoor = VMRI.mat*(Tmat\(VCT.mat\[calCTCoor ones(size(calCTCoor,1),1)]'));
    MRICoor = MRICoor(1:3,:)';

end


function allCalibCoord = centroidCalibration(allCoord,vol,vox2ras)
% Applying a barycenter calibration to the coord

% threshholding CT image
vol(vol<1000) = 0;

% iteration radius
cubic_r = 3;

% maximum adjustment range
max_adr = 2;
% maximum contacts interval limit
con_inter = [3 4];
% minimal iteration deviation limit
emin = 0.2;
% max iteration time
maxiter = 100;

allCalibCoord = zeros(size(allCoord));
for ic = 1:size(allCoord,1)
    
coord = allCoord(ic,:);    
    
calibCoord = [];
coord_vol = vox2ras\[coord ones(size(coord,1),1)]';
coord_vol = round(coord_vol(1:3,:)');

%% single input coord without neighbouring constraint
if size(coord,1) == 1
    
    initp_vol = coord_vol;
    
    initp_ras = coord';
    
    cm_vol = initp_vol;
    cm = initp_ras;
    cmo = initp_ras;
    
    % barycentric iteration
    iterator = 1;
    while iterator <= maxiter && cubic_r <= 6
        
        
        cubic_x = (cm_vol(1)-cubic_r:cm_vol(1)+cubic_r);
        cubic_y = (cm_vol(2)-cubic_r:cm_vol(2)+cubic_r);
        cubic_z = (cm_vol(3)-cubic_r:cm_vol(3)+cubic_r);
        
        Mx = 0;
        My = 0;
        Mz = 0;
        Mxyz = 0;
        for i = cubic_x
            for j = cubic_y
                for k = cubic_z
                    Mxyz = Mxyz + vol(i+1,j+1,k+1);
                    Mx = Mx + vox2ras(1,:)*[i;j;k;1]*vol(i+1,j+1,k+1);
                    My = My + vox2ras(2,:)*[i;j;k;1]*vol(i+1,j+1,k+1);
                    Mz = Mz + vox2ras(3,:)*[i;j;k;1]*vol(i+1,j+1,k+1);
                end
            end
        end
        
        % checking for reliability of centroid coordinate
        if all([Mx My Mz Mxyz])
            cm = [Mx/Mxyz;My/Mxyz;Mz/Mxyz];
        else
            cubic_r = cubic_r+1;
        end
        
        if norm(cm-cmo) >= emin && (iterator/(10*cubic_r)) <= cubic_r
            cmo = cm;
            cm_vol = round(vox2ras\[cm;1]);
        else if norm(cm-cmo) >= emin && (iterator/(10*cubic_r)) > cubic_r
                cmo = cm;
                cm_vol = round(vox2ras\[cm;1]);
                cubic_r = cubic_r+1;
            else break
            end
        end
        
        iterator = iterator+1;
    end
    
    calibCoord = cm';
    
end


%% double input coord with neighbouring constraint
if size(coord,1) == 2
    
    initp_vol = coord_vol(2,:);
    
    initp_ras = coord(2,:)';
    
    cm_vol = initp_vol;
    cm = initp_ras;
    cmo = initp_ras;
    
    % barycentric iteration
    iterator = 1;
    while iterator <= maxiter && cubic_r <= 2
        
        
        cubic_x = (cm_vol(1)-cubic_r:cm_vol(1)+cubic_r);
        cubic_y = (cm_vol(2)-cubic_r:cm_vol(2)+cubic_r);
        cubic_z = (cm_vol(3)-cubic_r:cm_vol(3)+cubic_r);
        
        Mx = 0;
        My = 0;
        Mz = 0;
        Mxyz = 0;
        for i = cubic_x
            for j = cubic_y
                for k = cubic_z
                    Mxyz = Mxyz + vol(i+1,j+1,k+1);
                    Mx = Mx + vox2ras(1,:)*[i;j;k;1]*vol(i+1,j+1,k+1);
                    My = My + vox2ras(2,:)*[i;j;k;1]*vol(i+1,j+1,k+1);
                    Mz = Mz + vox2ras(3,:)*[i;j;k;1]*vol(i+1,j+1,k+1);
                end
            end
        end
        
        % checking for reliability of centroid coordinate
        if all([Mx My Mz Mxyz])
            cm = [Mx/Mxyz;My/Mxyz;Mz/Mxyz];
        else break
        end
        
        % constraining the maximum adjustment range
        if norm(cm'-coord(1,:)) > con_inter(2) ...
                || norm(cm'-coord(1,:)) < con_inter(1) ...
                || norm(cm-initp_ras) > max_adr
            cm = cmo;
            break
        end
        
        if norm(cm-cmo) >= emin && (iterator/(10*cubic_r)) <= cubic_r
            cmo = cm;
            cm_vol = round(vox2ras\[cm;1]);
        else
            if norm(cm-cmo) >= emin && (iterator/(10*cubic_r)) > cubic_r
                cmo = cm;
                cm_vol = round(vox2ras\[cm;1]);
                cubic_r = cubic_r+1;
            else
                break
            end
        end
        
        iterator = iterator+1;
    end
    
    calibCoord = cm';
    
    
end

allCalibCoord(ic,:) = calibCoord;
end


