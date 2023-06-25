%% script for CW re-reference

% this script calculates the Closest White matter re-reference, based on the
% paper on Neuroimage  Arnulfo et al., <Phase and amplitude correlations in resting-state activity in human stereotactical EEG recordings>

% required input and corresponding size: rawsignal: channel * frames,
% elecposRas: channel * 3, badchannels : vector, fsSubPath: char, the
% .pial and .smoothwm file directory 
% outputs: CW_rereferenced_signal,  refContact: white matter reference
% contact for each contact, refDist: distance between contacts and
% their white matter reference contacts, gmpiPool: gmpi values for each
% contact, usage: becomes 0 if any contact has no appropriate white matter
% reference.


function [cwReferencedSignal,refContact,refDist,usage,gmpiPool] = CWReref(rawsignal,elecposRas,badchannels,fsSubPath)


wm_thresh = -0.7; % calculate from our image data of 21 patients (-0.7411 actual)

usage=1;


 % NOTE: make sure the vertex_coords is in tkr space!!
 
[vertex_coords1, faces1]=read_surf([fsSubPath '/surf/lh.pial']);
faces1 = faces1 + 1;

[vertex_coords2, faces2]=read_surf([fsSubPath '/surf/rh.pial']);
faces2 = faces2 + 1;

faces_p = [faces1;faces2+size(vertex_coords1,1)];

vertex_coords_p = [vertex_coords1;vertex_coords2];

[vertex_coords1, faces1]=read_surf([fsSubPath '/surf/lh.smoothwm']);
faces1 = faces1 + 1;

[vertex_coords2, faces2]=read_surf([fsSubPath '/surf/rh.smoothwm']);
faces2 = faces2 + 1;

faces_w = [faces1;faces2+size(vertex_coords1,1)];

vertex_coords_w = [vertex_coords1;vertex_coords2];






%% calculate GMPI value for all contacts
elecposRas = elecposRas(1:size(rawsignal,1),:);
gmpiPool = ones(1,size(elecposRas,1));


for contact_ind = 1:size(elecposRas,1)
    
current_elec_coor = elecposRas(contact_ind,:);

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

end



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
       
       % what if it is a bad channel??
       
      
end


cwReferencedSignal = rawsignal - rawsignal(refContact,:);


refDist = PairDistance(elecposRas,elecposRas(refContact,:)); 

cwReferencedSignal(badchannels,:) = 0;

refDist(badchannels)=[];
    

end


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






