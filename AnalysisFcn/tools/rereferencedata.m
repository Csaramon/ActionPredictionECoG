function  varargout = rereferencedata(ddata,re_reference,val_elecs,NumArray,ref_parameter,surfdir)
% Reference data using different method
%

FiltReref_data = zeros(size(ddata));
distances_dist =[];
switch re_reference
    
    case 'no_reference'
        FiltReref_data  = ddata ;
        
    case 'common_average'
        % common average reference,if this is depth electrode, no need to re-reference
        CAR = mean(ddata(val_elecs,:));  % common average,use all the valid electrode
        for ielec = val_elecs
            FiltReref_data(ielec,:) = ddata(ielec,:) - CAR;  %=%  common average reference
        end
        
    case 'strip_average'
        % average reference among each electrode strip
        StrNum = cumsum(NumArray);
        StrNum = [0 StrNum];
        for s = 1:(length(StrNum)-1)
            A = [StrNum(s)+1:StrNum(s+1)];
            StrVal{s} = A(ismember(A,val_elecs));  % valid contacts in each strip
        end
                
        for p = 1:length(StrVal)
            ee = StrVal{p};
            SAR = mean(ddata(ee,:));
            FiltReref_data(ee,:) = ddata(ee,:) - repmat(SAR,[length(ee) 1]);
        end
        
    case 'surface_depth_ECoG'
        % average reference among electrode with strips and patches
            depth_strip = ref_parameter;
        StrNum = cumsum(NumArray);
        StrNum = [0 StrNum];
        s = 1:(length(StrNum)-1);
        surface_patch = setdiff(s,depth_strip);
        
        % common reference for surface electrodes
        surface_contact = [];
        for ss = surface_patch
            surface_contact = [surface_contact StrNum(ss)+1:StrNum(ss+1)];  %% all the surface contacts
        end
        surface_valid = surface_contact(ismember(surface_contact,val_elecs)); %% valid surface contacts
        surface_ave = mean(ddata(surface_valid,:));  %% average of the valid surface contacts
        
        for e = surface_valid
            FiltReref_data(e,:) = ddata(e,:) - surface_ave;
        end
        
        % strip reference for depth electrodes
        for s = depth_strip
            A = [StrNum(s)+1:StrNum(s+1)];
            StrVal{s} = A(ismember(A,val_elecs));  % valid contacts in this strip
        end
        
        for p = 1:length(StrVal)
            ee = StrVal{p};
            SAR = mean(ddata(ee,:));
            FiltReref_data(ee,:) = ddata(ee,:) - repmat(SAR,[length(ee) 1]);
        end
        
    case  'WhiteMatter_ComAve'
        % only common average of the white matter contacts as reference
        WhitMattCont = ref_parameter;
        wmvalid = ismember(WhitMattCont,val_elecs);
        refcont = WhitMattCont(wmvalid);
        
        CAR = mean(ddata(refcont,:));  % common average,use all the valid electrode
        for ielec = val_elecs
            FiltReref_data(ielec,:) = ddata(ielec,:) - CAR;  %  common reference
        end
        
        
    case  'WhiteMatter_StrAve'
        % only strip average of the white matter contacts as reference
        WhitMattCont = ref_parameter;
        wmvalid = ismember(WhitMattCont,val_elecs);
        refcont = WhitMattCont(wmvalid);
        StrNum = cumsum(NumArray);
        StrNum = [0 StrNum];
        for s = 1:(length(StrNum)-1)
            A = [StrNum(s)+1:StrNum(s+1)];
            StrVal{s} = A(ismember(A,refcont));  % valid contacts in this strip
        end
        
        for p = 1:length(StrVal)
            ee = StrVal{p};
            SAR = mean(ddata(ee,:));
            FiltReref_data(ee,:) = ddata(ee,:) - repmat(SAR,[length(ee) 1]);
        end
        
        
    case 'Bipolar'
        % transform to bipolar configration
        % bipolar every electrode would lose a contact which might cause trouble
        % for later analysis. To offset this, the last bipolar contact is
        % replicated and thus make number of contacts on every electrode unchanged.
        sumelecnum = [0 cumsum(NumArray)];
        bipolardata = [];
        for i =1:length(NumArray)
            
            jianshu_vector = ddata([sumelecnum(i)+1: (sumelecnum(i+1)-1) sumelecnum(i+1)-1],: );
            
            beijianshu_vector = ddata([sumelecnum(i)+2:sumelecnum(i+1) sumelecnum(i+1)],: );
            
            bipolar_this_elec = beijianshu_vector - jianshu_vector;
            
            bipolardata = [bipolardata;bipolar_this_elec];
            
        end
        FiltReref_data = bipolardata;
        
        
    case 'Closest WhiteMatter'
        
        badchannels = setdiff(1:sum(NumArray),val_elecs);
        [FiltReref_data,reference_pool,distances_dist,~,~] = cw_calculation(ddata,ref_parameter,badchannels,surfdir);
        
end

varargout{1} = FiltReref_data;
varargout{2} = distances_dist;
end
