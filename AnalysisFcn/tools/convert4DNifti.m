niiFileName = '~/Desktop/CingulateOrbitoFrontal_thr25_1mm.nii.gz';


aparc_nii = load_nifti(niiFileName);
aparc_vol = aparc_nii.vol;
newVol = zeros(size(aparc_vol,1),size(aparc_vol,2),size(aparc_vol,3));


for i = 1:size(newVol,1)
    for j = 1:size(newVol,2)
        for k = 1:size(newVol,3)
            
            [Vmax,Imax] = max(squeeze(aparc_vol(i,j,k,:)));
            if Vmax > 0
                newVol(i,j,k) = Imax;

            end
        
        end
    end
end

aparc_nii.vol = newVol;

aparc_nii = save_nifti(aparc_nii,'~/Desktop/PFCMaxProb.nii');