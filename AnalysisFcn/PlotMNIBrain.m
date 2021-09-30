%% plot MNI152 pial surface
vox2ras = [-1,0,0,127;0,0,1,-145;0,-1,0,147;0,0,0,1];
ras_tkr2vox = [-1,0,0,128;0,0,-1,128;0,1,0,128;0,0,0,1];
[~,MNI152path] = unix('. ~/.zshrc;echo $SUBJECTS_DIR ');
MNI152path = [MNI152path(1:end-1) filesep 'cvs_avg35_inMNI152'];
% set figure handle
stereoview = figure('Name','3D View','NumberTitle','of');
stereoaxes = axes('Parent',stereoview);
hold(stereoaxes,'on')
axis(stereoaxes,'equal');
axis(stereoaxes,'off');
% plot pial surface
[vertex_coords1, faces1]=read_surf([MNI152path '/surf/lh.pial']);
faces1 = faces1 + 1;
[vertex_coords2, faces2]=read_surf([MNI152path '/surf/rh.pial']);
faces2 = faces2 + 1;
% stereonum = [size(vertex_coords1,1) size(faces1,1)];
% faces = [faces1;faces2+size(vertex_coords1,1)];
% vertices = [vertex_coords1;vertex_coords2];
% left only
faces = [faces1];
vertices = [vertex_coords1];

vertices = vox2ras*ras_tkr2vox*[vertices ones(size(vertices,1),1)]';
vertices = vertices(1:3,:)';
stereo = patch(struct(...
    'vertices', vertices, 'faces', faces), ...
    'Parent',stereoaxes, ...
    'FaceColor',[230,228,216]./255, ...
    'FaceAlpha',0.25, ...
    'EdgeColor', 'none');

brainLight = camlight;
view(-90,0)
set(brainLight,'position',get(gca,'cameraPosition'))
lighting gouraud;
material dull

