HCP    ='/home/essicd/storage/data/HCP/Unrelated_80';
Outdir ='/storage/maullz/Contour_Inference_2018/HCP_analysis_data_test';

if ~isdir(Outdir)
    mkdir(Outdir)
end

% Cope files
Copes     = dir([HCP '/*/WM/Level2/cope11.feat/cope1*']);
Varcopes  = dir([HCP '/*/WM/Level2/cope11.feat/varcope1*']);
Tstats    = dir([HCP '/*/WM/Level2/tstat11.nii.gz']);

% Copying to Outdir
for i=1:length(Copes)
    dir_split = strsplit(Copes(i).folder,'/');
    sub_id    = dir_split{7}; 
    copyfile(fullfile(Copes(i).folder,Copes(i).name), fullfile(Outdir, ['cope_' sub_id '.nii.gz'])); 
    command   = ['fslmaths ' fullfile(Varcopes(i).folder,Varcopes(i).name) ' -bin ' fullfile(Outdir, ['mask_' sub_id '.nii.gz'])];
    system(command);    
    % Adding additional smoothing to the Cope image so that final image has 6mm FWHM smoothing. Using the squation FWHM = sigma*sqrt(8ln2) with FWHM = sqrt(6^2 - 4^2) 
    command   = ['fslmaths ' fullfile(Outdir, ['cope_' sub_id '.nii.gz']) ' -kernel gauss 1.9 -fmean ' fullfile(Outdir, ['smooth_cope_' sub_id '.nii.gz'])];
    system(command);
end

% Concatenate copes to 4D file
out_copes = cellstr(spm_select('FPList', Outdir, 'cope.*\.nii.gz$'));
command = strcat('fslmerge -t', {' '}, fullfile(Outdir,'copes.nii.gz'), {' '}, strjoin(out_copes'));
system(char(command));

% Concatenate smooth copes to 4D file
out_copes = cellstr(spm_select('FPList', Outdir, 'smooth_cope.*\.nii.gz$'));
command = strcat('fslmerge -t', {' '}, fullfile(Outdir,'smooth_copes.nii.gz'), {' '}, strjoin(out_copes'));
system(char(command));

% Create group mask from the original, unsmoothed copes
out_masks = cellstr(spm_select('FPList', Outdir, 'mask.*\.nii.gz$'));
command = strcat('fslmerge -t', {' '}, fullfile(Outdir,'group_mask.nii.gz'), {' '}, strjoin(out_masks'));
system(char(command));

command = ['fslmaths ' fullfile(Outdir,'group_mask.nii.gz') ' -Tmin ' fullfile(Outdir,'group_mask.nii.gz')];
system(command); 

% Unzipping the copes so the data can be analyzed in SPM
for i=1:length(Copes)
    dir_split = strsplit(Copes(i).folder,'/');
    sub_id    = dir_split{7}; 
    command   = ['gunzip ' fullfile(Outdir, ['smooth_cope_' sub_id '.nii.gz'])];
    system(command);
end 

 