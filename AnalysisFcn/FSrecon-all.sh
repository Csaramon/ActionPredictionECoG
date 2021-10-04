# script for free surfer reconstruction

export SUBJECTS_DIR=/data00/Chaoyi/ActionPrediction/Patient13/

mri_convert -c -oc 0 0 0 $SUBJECTS_DIR/T1.nii $SUBJECTS_DIR/T1.nii

recon-all -s FsRecon -i $SUBJECTS_DIR/T1.nii  -all

mri_cvs_register --mov FsRecon --template cvs_avg35_inMNI152 --templatedir /usr/local/freesurfer/subjects/
