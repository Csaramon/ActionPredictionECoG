function CovertCoordinates
[CTfilename, CTpathname, ~] = uigetfile(['/Users/qinchaoyi/Desktop/*.nii'],'Choose CT scan to register');
[MRIfilename, MRIpathname, ~] = uigetfile(['/Users/qinchaoyi/Desktop/*.nii'],'Choose template MRI scan');
[ROIpathname] = uigetdir('/Users/qinchaoyi/Desktop', 'Pick the directory containing the ROI labels');

VCT = spm_vol([CTpathname CTfilename]);
CTVOL = spm_read_vols(VCT);
VMRI = spm_vol([MRIpathname MRIfilename]);

ROIFiles = dir([ROIpathname filesep '*.label']);

F = spm_figure('Create','Graphics','coregistration',1);
flags.graphics = 1;
x = spm_coreg(VMRI,VCT,flags);
Tmat = VCT.mat\spm_matrix(x(:)')*VMRI.mat;

allMRICoor = [];
for iroi = 1:numel(ROIFiles)
    
    CTCoor = read_label([],[ROIFiles(iroi).folder filesep ROIFiles(iroi).name]);
    CTCoor = CTCoor(:,2:4);
    
    calCTCoor = centroidCalibration(CTCoor,CTVOL,VCT.mat);
    
    
    MRICoor = VMRI.mat*(Tmat\(VCT.mat\[calCTCoor ones(size(calCTCoor,1),1)]'));
    MRICoor = MRICoor(1:3,:)';
    
    allMRICoor{iroi,1} = MRICoor;
end

uiwait(F)

%% determine the electrodes' order

[datafilename, datapathname, ~] = uigetfile(['/Users/qinchaoyi/Desktop/*.edf'],'Choose edf to get hdr info');

hdr = edfread([datapathname datafilename]);

newLabel = {};
for ich = 1:numel(hdr.label)
    nameInd = find(hdr.label{ich}=='-');
    newLabel{ich} = hdr.label{ich}(1:nameInd-1);
end

uniLabel = unique(newLabel,'stable')';
tchan = tabulate(newLabel);

global gdata

gdata.uniLabel = uniLabel;
gdata.allMRICoor = allMRICoor;
gdata.newuniLabel = uniLabel(1:numel(ROIFiles));

gdata.h = figure('Name','All clusters','NumberTitle','off','color',[1 1 1]);
gdata.ha =axes('Parent',gdata.h,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'ztick',[],'zticklabel',[]);
axis equal
axis off
hold all
rotate3d on

% configuring the ui figure
gdata.alpha_text = uicontrol(gdata.h,'style','text', ...
    'backgroundcolor',[1 1 1], ...
    'units', 'normalized', ...
    'position',[0.92 0.4 0.08 0.04], ...
    'string','Surface');

gdata.alpha_value = uicontrol(gdata.h,'style','edit', ...
    'backgroundcolor',[1 1 1], ...
    'units', 'normalized', ...
    'position',[0.93 0.36 0.06 0.04], ...
    'string','0.15', ...
    'callback','global gdata;set(gdata.stereo,''facealpha'',str2double(get(gdata.alpha_value,''string'')));set(gdata.light,''position'',get(gca,''cameraPosition''))');

gdata.panel = uipanel(gdata.h,'Title','Choosing electrode clusters','FontSize',12,...
    'BackgroundColor',[242 255 198]/255,...
    'Position',[0 0 1 0.2]);

gdata.continue = uicontrol(gdata.h,'style','pushbutton', ...
    'units', 'normalized', ...
    'String','Continue', ...
    'position',[0.92 0.28 0.08 0.06], ...
    'callback','uiresume(gcbf)');

% show fs surface
[vertex_coords1, faces1]=read_surf([MRIpathname(1:end-4) '/surf/lh.pial']);
faces1 = faces1 + 1;
[vertex_coords2, faces2]=read_surf([MRIpathname(1:end-4) '/surf/rh.pial']);
faces2 = faces2 + 1;
gdata.faces = [faces1;faces2+size(vertex_coords1,1)];
gdata.vertices = [vertex_coords1;vertex_coords2];
gdata.stereo = patch(struct(...
    'vertices', gdata.vertices, 'faces', gdata.faces), ...
    'Parent',gdata.ha, ...
    'FaceColor',[230,228,216]./255, ...
    'FaceAlpha',0.15, ...
    'EdgeColor', 'none');
material dull
gdata.light = camlight;
set(gdata.light,'position',get(gca,'cameraPosition'))
lighting gouraud;


ypos = 1;
yrow = ceil(numel(ROIFiles)/4);
if ~isfield(gdata,'clusterpoints') || isempty(gdata.clusterpoints)
    gdata.clusterpoints = cell(1,numel(ROIFiles));
end
subColor = distinguishable_colors(numel(ROIFiles));
for iarray = 1:numel(ROIFiles)
    
    xpos = mod(iarray,4);
    if xpos == 0
        xpos = 4;
    end
    
    tempCoor = allMRICoor{iarray};
    for ie = 1:size(tempCoor,1)
        
        scatter3(gdata.ha,tempCoor(ie,1),tempCoor(ie,2),tempCoor(ie,3),60,subColor(iarray,:),'fill');
    end
    text(tempCoor(end,1),tempCoor(end,2)+2,tempCoor(end,3)+1,num2str(iarray),'Fontsize',25,'Color',subColor(iarray,:))
    
    ecluster_text(iarray) = uicontrol(gdata.panel,'style','text', ...
        'string', num2str(iarray),...
        'units', 'normalized', ...
        'FontSize',12, ...
        'BackgroundColor',[242 255 198]/255,...
        'position',[0.01+0.25*(xpos-1) 0.98-0.96/yrow*ypos 0.04 0.96/yrow]);
    
    gdata.ecluster(iarray) = uicontrol(gdata.panel,'style','popupmenu', ...
        'UserData',iarray, ...
        'Value',iarray,...
        'string', uniLabel,...
        'backgroundcolor',[1 1 1], ...
        'units', 'normalized', ...
        'position',[0.05+0.25*(xpos-1) 0.99-0.96/yrow*ypos 0.2 0.96/yrow], ...
        'callback',@ChooseElectrode);
    
    
    if xpos == 4
        ypos = ypos+1;
    end
    
end

% continue the function
uiwait(gdata.h)
try
    close(gdata.h)
end

tempLabel = [];
tempCoor = [];
for ia = 1:numel(gdata.newuniLabel)

    il = strncmp(hdr.label,gdata.newuniLabel{ia},length(gdata.newuniLabel{ia}));
    tempLabel = [tempLabel;hdr.label(il)'];
    if sum(il) < size(gdata.allMRICoor{ia},1)
        warndlg(['Number of contacts are mismatch (' num2str(size(gdata.allMRICoor{ia},1)) ' coordinates compare to ',...
            num2str(sum(il)) ' labels in ' gdata.newuniLabel{ia} '),'  ...
            'outermost contact coordinates will be removed!'])
    elseif sum(il) > size(gdata.allMRICoor{ia},1)
        errordlg(['Number of contacts are mismatch (' num2str(size(gdata.allMRICoor{ia},1)) ' coordinates compare to ',...
            num2str(sum(il)) ' labels in ' gdata.newuniLabel{ia} '),'  ...
            'Please check the coordinates on the CT!'])
        return
    end
    tempCoor = [tempCoor;gdata.allMRICoor{ia}(1:sum(il),:)];
end

% sort all labels and cooridnates to their original orders
newLabel = [];
newCoor = [];
for il = 1:numel(hdr.label)
    tempInd = strcmp(hdr.label(il),tempLabel);
 if any(tempInd)
    newLabel = [newLabel;tempLabel(tempInd)];
    newCoor = [newCoor;tempCoor(tempInd,:)];
 end
end
%% transforma coordiantes to MNI space

% generate temporary coordinates file
fid=fopen([ROIpathname filesep 'tempcoordiantes.txt'],'w+');
for i = 1:size(newCoor,1)
    fprintf(fid,'%f\t%f\t%f\t%f\n',newCoor(i,:),i);
end
fclose(fid);

% generate roi file to warp in MNI152 space
setbash = ['. ~/.zshrc;PATH=$PATH:/Users/qinchaoyi/Documents/MATLAB/toolbox/FIELD/addon/;'];
gen_roi2warp_code = [setbash ...
    '3dUndump -master ' MRIpathname '/T1.nii -orient LPI ' ...
    '-srad 3 -xyz -prefix ' ROIpathname filesep '/elec_rois_to_warp_in_mni152.nii ' ...
    ROIpathname filesep 'tempcoordiantes.txt'];
unix(gen_roi2warp_code)

% generate roi in MNI152 space
[~,MNI152path] = unix('. ~/.zshrc;echo $SUBJECTS_DIR ');
MNI152path = [MNI152path(1:end-1) filesep 'cvs_avg35_inMNI152'];
genMNI152code = [setbash 'mri_vol2vol --targ ' MNI152path '/mri/T1.mgz --m3z '  ...
    MRIpathname(1:end-4) '/cvs/final_CVSmorph_tocvs_avg35_inMNI152.m3z --noDefM3zPath --mov '  ...
    ROIpathname '/elec_rois_to_warp_in_mni152.nii --o ' ...
    ROIpathname '/elec_rois_MNI152.nii --interp nearest --no-save-reg'];

unix(genMNI152code)

mni_info = load_nifti([ROIpathname '/elec_rois_MNI152.nii']);
mni_vol = mni_info.vol;

volCoordinate = zeros(size(newCoor));
for ii = 1:size(newCoor,1)
    ind = find(mni_vol==ii);
    [sub1,sub2,sub3] = ind2sub(size(mni_vol),ind);
    volCoordinate(ii,:) = [mean(sub1) mean(sub2) mean(sub3)];
end
MNICoordinate = mni_info.vox2ras*[volCoordinate ones(size(volCoordinate,1),1)]';
MNICoordinate = MNICoordinate(1:3,:)';
save([ROIpathname '/elec.mat'],'MNICoordinate','newCoor','newLabel')
clear global
end


%% subfunctions


function ChooseElectrode(hObject,~)
global gdata

gdata.newuniLabel{hObject.UserData} = hObject.String{hObject.Value};
% gdata.newMRICoor{hObject.UserData} = gdata.allMRICoor{hObject.Value};
% if size(gdata.allMRICoor{hObject.Value},1) > gdata.numContact{hObject.Value}
%     warndlg(['Number of contacts are mismatch (' num2str(size(gdata.allMRICoor{hObject.Value},1)) ' coordinates compare to ',...
%         num2str(gdata.numContact{hObject.Value}) ' labels),'  ...
%         'outermost contact coordinates will be removed!'])
%     gdata.newMRICoor{hObject.UserData} = gdata.allMRICoor{hObject.UserData}(1:gdata.numContact{hObject.Value},:);
% end
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

end


