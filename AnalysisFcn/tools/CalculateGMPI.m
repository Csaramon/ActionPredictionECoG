function CalculateGMPI

% base path to subjects' folder
base_path = '/Users/apple/Desktop/BistableProject/MyAnalysis/';
% specify the subjects to use (subject folder arranged in freesurfer structure)
allsub = {'chenqiuhong','houhaoran','liaojianlin','wangqiangqiang','zhaojie', ...
    'chenhuaping','chensa','dongfenglian','fanjiang','guobin',...
    'langweichuang','liqiang','liuzhiming','luping','peijian','tianyongqiang', ...
    'wangjunjie','wangjunlu','wangmingyue','zhangfeng','zhaozhichao'}; % all sub


% section1: Calculate GMPI threshold
pWM = 0.05; % possibility of electrode in WM using calculated GMPI threshold

gmpiAll = [];
aparcIndexAll = [];
for isub = 1:numel(allsub)
    fprintf(['Currently processing: sub ' num2str(isub) ' in ' num2str(numel(allsub))  ' id: ' allsub{isub} '\n']);
    subpath = [base_path allsub{isub}];
    surfDir = [subpath '/surf/'];
    
    % load coordinates
    aa = load([subpath '/brain3D/autocoordinates.mat']);
    cc = fieldnames(aa);
    elecposRaw = aa.(cc{1});
    elecposRas = elecposRaw(:,end-2:end);
    aparcNii = load_nifti([subpath '/mri/aparc+aseg.nii']);
    aparcVol = int32(aparcNii.vol);
    
    vox2ras = inv(aparcNii.vox2ras);
    ras2vox = inv(aparcNii.vox2ras);
    vox2ras_tkr = [ -1 0 0 128; 0 0 1 -128; 0 -1 0 128; 0 0 0 1];
    
    elecposVol = ras2vox*[elecposRas(:,1:3) ones(size(elecposRas,1),1)]';
    elecposVol = elecposVol(1:3,:)';
    elecposTkr = vox2ras_tkr*ras2vox*[elecposRas ones(size(elecposRas,1),1)]';
    elecposTkr = elecposTkr(1:3,:)';
    % load surface files
    % NOTE: make sure the vertex_coords is in tkr space!!
    [vertex_coords1, faces1]=read_surf([surfDir 'lh.pial']);
    faces1 = faces1 + 1;
    
    [vertex_coords2, faces2]=read_surf([surfDir 'rh.pial']);
    faces2 = faces2 + 1;
    
    faces_p = [faces1;faces2+size(vertex_coords1,1)];
    
    vertex_coords_p = [vertex_coords1;vertex_coords2];
    
    [vertex_coords1, faces1]=read_surf([surfDir 'lh.smoothwm']);
    faces1 = faces1 + 1;
    
    [vertex_coords2, faces2]=read_surf([surfDir 'rh.smoothwm']);
    faces2 = faces2 + 1;
    
    faces_w = [faces1;faces2+size(vertex_coords1,1)];
    
    vertex_coords_w = [vertex_coords1;vertex_coords2];
    
    %% calculate GMPI value for all contacts
    gmpiPool = zeros(1,size(elecposTkr,1));
    aparcIndex = zeros(1,size(elecposTkr,1));
    
    
    for contact_ind = 1:size(elecposTkr,1)
        
        current_elec_coor = elecposTkr(contact_ind,:);
        
        dist_vec_p = one_to_all_dist(current_elec_coor,vertex_coords_p);
        
        dist_vec_w = one_to_all_dist(current_elec_coor,vertex_coords_w);
        
        [~,nearest_p] = min(dist_vec_p);
        
        [~,nearest_w] = min(dist_vec_w);
        
        p_coor = vertex_coords_p(nearest_p,:);
        
        w_coor = vertex_coords_w(nearest_w,:);
        
        wc = current_elec_coor - w_coor;
        
        wp = p_coor - w_coor;
        
        norm_wp = norm(wp);
        
        gmpi = (wc*wp')./norm_wp;
        
        gmpiPool(contact_ind)=gmpi;
        
        
        ccVol = elecposVol(contact_ind,:);
        
        aparcIndex(contact_ind) = aparcVol(round(ccVol(1)),round(ccVol(2)),round(ccVol(3)));
        
    end
    
    gmpiAll = [gmpiAll,gmpiPool];
    aparcIndexAll = [aparcIndexAll,aparcIndex];

end

% plot GM,WM,Unknown GMPI histogram
figure;
hall = histogram(gmpiAll);
hold on

gmpiWM = gmpiAll(aparcIndexAll==2 | aparcIndexAll==41);
gmpiUN = gmpiAll(aparcIndexAll==0);
gmpiGM = gmpiAll(aparcIndexAll~=2 & aparcIndexAll~=41 & aparcIndexAll~=0);

gmpiGMsort = sort(gmpiGM);
gmpiGM005 = gmpiGMsort(round(numel(gmpiGMsort)*pWM));
hwm = histogram(gmpiWM);
hun = histogram(gmpiUN);
hgm = histogram(gmpiGM);
h005 = plot([gmpiGM005 gmpiGM005],[0 max(hall.Values)],'LineWidth',2,'Color',[0 0 0]);
legend([hall,hwm,hgm,hun,h005],['ALL'],['WM'],['GM'],['Unknown'],['GMPI = ' num2str(gmpiGM005) '(p=' num2str(pWM) ')'])
xlim([-8 8])
xlabel('GMPI')
ylabel('Count')

fprintf(['GMPI threshold for the current dataset is: ' num2str(gmpiGM005) '\n'])



%% section2: calculate and plot correlation between GMPI threshold and elec2ref
% distance
wmthreshAll = linspace(-2,1,20);
refDistYmean = zeros(size(wmthreshAll));
refDistYstd = zeros(size(wmthreshAll));
for it = 1:numel(wmthreshAll)
    fprintf(['Currently processing datapoints: ' num2str(it) ' in ' num2str(numel(wmthreshAll)) '\n']);
    wm_thresh = wmthreshAll(it);
    
    
    for isub = 1:numel(allsub)
        
        subpath = [base_path allsub{isub}];
        surf_dir = [subpath '/surf/'];
        
        % load coordinates
        aa = load([subpath '/brain3D/autocoordinates.mat']);
        cc = fieldnames(aa);
        elecposRaw = aa.(cc{1});
        elecposRas = elecposRaw(:,end-2:end);
        aparcNii = load_nifti([subpath '/mri/aparc+aseg.nii']);
        aparcVol = int32(aparcNii.vol);
        
        vox2ras = inv(aparcNii.vox2ras);
        ras2vox = inv(aparcNii.vox2ras);
        vox2ras_tkr = [ -1 0 0 128; 0 0 1 -128; 0 -1 0 128; 0 0 0 1];
        
        elecposVol = ras2vox*[elecposRas(:,1:3) ones(size(elecposRas,1),1)]';
        elecposVol = elecposVol(1:3,:)';
        elecposTkr = vox2ras_tkr*ras2vox*[elecposRas ones(size(elecposRas,1),1)]';
        elecposTkr = elecposTkr(1:3,:)';
        % load surface files
        % NOTE: make sure the vertex_coords is in tkr space!!
        [vertex_coords1, faces1]=read_surf([surf_dir 'lh.pial']);
        faces1 = faces1 + 1;
        
        [vertex_coords2, faces2]=read_surf([surf_dir 'rh.pial']);
        faces2 = faces2 + 1;
        
        faces_p = [faces1;faces2+size(vertex_coords1,1)];
        
        vertex_coords_p = [vertex_coords1;vertex_coords2];
        
        [vertex_coords1, faces1]=read_surf([surf_dir 'lh.smoothwm']);
        faces1 = faces1 + 1;
        
        [vertex_coords2, faces2]=read_surf([surf_dir 'rh.smoothwm']);
        faces2 = faces2 + 1;
        
        faces_w = [faces1;faces2+size(vertex_coords1,1)];
        
        vertex_coords_w = [vertex_coords1;vertex_coords2];
        
        %% calculate GMPI value for all contacts
        gmpiPool = zeros(1,size(elecposTkr,1));
        aparcIndex = zeros(1,size(elecposTkr,1));
        
        
        for contact_ind = 1:size(elecposTkr,1)
            
            current_elec_coor = elecposTkr(contact_ind,:);
            
            dist_vec_p = one_to_all_dist(current_elec_coor,vertex_coords_p);
            
            dist_vec_w = one_to_all_dist(current_elec_coor,vertex_coords_w);
            
            [~,nearest_p] = min(dist_vec_p);
            
            [~,nearest_w] = min(dist_vec_w);
            
            p_coor = vertex_coords_p(nearest_p,:);
            
            w_coor = vertex_coords_w(nearest_w,:);
            
            wc = current_elec_coor - w_coor;
            
            wp = p_coor - w_coor;
            
            norm_wp = norm(wp);
            
            gmpi = (wc*wp')./norm_wp;
            
            gmpiPool(contact_ind)=gmpi;
            
            
            ccVol = elecposVol(contact_ind,:);
            
            aparcIndex(contact_ind) = aparcVol(round(ccVol(1)),round(ccVol(2)),round(ccVol(3)));
            
        end
        
        gmpiAll = [gmpiAll,gmpiPool];
        aparcIndexAll = [aparcIndexAll,aparcIndex];
        
        %% rereference to CW scheme
        
        [distance_vector,mean_dist, std_dist, distance_matrix] = point_dist(elecposRas);
        
        refContact = zeros(1,size(elecposRas,1));
        
        
        for contact_ind = 1:size(elecposRas,1)
            
            current_elec_coor = elecposRas(contact_ind,:);
            
            
            [~,ind] = sort(distance_matrix(contact_ind,:));
            
            ini_ind = 2;
            
            
            while gmpiPool(ind(ini_ind)) > wm_thresh || ismember(ind(ini_ind),badchannels)
                
                ini_ind = ini_ind +1;
                
                if ini_ind>=length(ind)
                    
                    usage = 0;
                    
                    break;
                end
                
            end
            
            reference_contact_ind = ind(ini_ind);
            
            refContact(contact_ind) = reference_contact_ind;
            
        end
        
        refDist = PairDistance(elecposRas,elecposRas(refContact,:));
        
        refDist(badchannels)=[];
        
        refDistAll = [refDistAll,refDist];
        
    end
    
    
    
    % plot GM,WM,Unknown GMPI distance correlation
    refDistWM = refDistAll(aparcIndexAll==2 | aparcIndexAll==41);
    refDistUN = refDistAll(aparcIndexAll==0);
    refDistGM = refDistAll(aparcIndexAll~=2 & aparcIndexAll~=41 & aparcIndexAll~=0);
    
    refDistYmean(it) = mean(refDistAll);
    refDistYstd(it) = std(refDistAll);
end

figure;
plot(wmthreshAll,refDistYmean)
hold on
plot(wmthreshAll,refDistYmean+refDistYstd)
plot(wmthreshAll,refDistYmean-refDistYstd)
h005 = plot([gmpiGM005 gmpiGM005],[0 35],'LineWidth',2,'Color',[0 0 0]);
xlabel('GMPI')
ylabel('Distance(mm)')


end % function end








% sub fucntions 
function dist_vec = PairDistance(lista,listb)

dist_vec = ones(1,size(lista,1));

for k = 1:size(lista,1)
    
    
    dist_vec(k) = norm(lista(k,:)-listb(k,:));
    
    
end


end


% calculate the normed distance between points represented by each line of an matrix.


function [distance_vector,mean_dist, std_dist, distance_matrix]= point_dist(m)

distance_matrix=zeros(size(m,1),size(m,1));

for ind_a = 1:size(m,1)
    
    for ind_b = 1:size(m,1)
        
        distance_matrix(ind_a,ind_b) = norm(m(ind_a,:) - m(ind_b,:));
        
        
    end
    
end


N = length(distance_matrix);

a = zeros([1,N*(N-1)/2]);

for i=1:N;
    l = (i-1)*(2*N-i)/2+1;
    a(l:(l+N-i-1)) = distance_matrix((i+1):end,i);
end


distance_vector = a;

mean_dist = mean(a);

std_dist = std(a);


end



% calculates the euclidean distance vector of a point to a point list
% size of list should be N * 3

function dist_vec = one_to_all_dist(point,list)

dist_mat = repmat(point,[size(list,1) 1])-list;

dist_vec = ones(1,size(list,1)).*10;


for k = 1:size(list,1)
    
    
    dist_vec(k) = norm(dist_mat(k,:));
    
    
end

end






