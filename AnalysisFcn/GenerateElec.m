
ft_defaults
clear;clc

subname = 'Patient13';
subpath = ['/Users/qinchaoyi/Desktop/ActionPrediction/' subname '/'];

%% append electrode coordinates infomation
% load individual coordinates
aa = load([subpath '/FsRecon/brain3D/projectcoordinates.mat']);
cc = fieldnames(aa);
elecposRaw = aa.(cc{1});

% load MNI152 coordinates
elecposMNI = importdata([subpath '/FsRecon/brain3D/MNI152_coordinates_ras.txt']);
elecposMNI = elecposMNI(:,1:3);

% load session data for elec labels
bb = load([subpath '/RawData/' subname '_SequentialActionMovieTask_A1.mat']);
dd = fieldnames(bb);
hdr = bb.(dd{1});

elec = [];
chanarray = [];
elec.elecpos = elecposRaw(:,end-2:end); % electrode coordinates in RAS space
elec.label   = hdr.label; % electrode labels
elec.chanpos = elec.elecpos; % channel coordinates in RAS space(for seeg the same as elecpos?)
elec.coordsys = 'acpc'; % the coordinate system
arraynum = unique(round(elecposRaw(:,2)/100)); % electrode number
for ia = arraynum'
    chanarray{ia} = find(round(elecposRaw(:,2)/100)==ia); % number of contacts among each
end
elec = ft_datatype_sens(elec); % add electrode infomation to the fieldtrip data structure
elec.chanarray = chanarray;
elec.elecposMNI = elecposMNI;



%% Generate electrode index


c1 = 'A';
cn1 = 52;

c2 = 'B';
cn2 = 60;

newlabel = [];
for i1 = 1:cn1
    newlabel = [newlabel;{[c1,num2str(i1)]}];
end

for i2 = 1:cn2
    newlabel = [newlabel;{[c2,num2str(i2)]}];
end

elec.label = newlabel;
elec.chanunit = hdr.chanunit(1:numel(newlabel));

%% save electrode coordinates infomation
uisave('elec', [subpath '/FsRecon/brain3D/ft_elec.mat']);

