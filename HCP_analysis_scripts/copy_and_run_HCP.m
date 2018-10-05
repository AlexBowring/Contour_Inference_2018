HCP    ='/home/essicd/storage/data/HCP/Unrelated_80';
Outdir ='/storage/maullz/Contour_Inference_2018/HCP_analysis_results';

if ~isdir(Outdir)
    mkdir(Outdir)
end

% Cope files
Copes     = dir([HCP '/*/WM/Level2/cope11.feat/cope1*']);
Varcopes  = dir([HCP '/*/WM/Level2/cope11.feat/varcope1*']);

% Copying to Outdir
for i=1:length(Copes)
    dir_split = strsplit(Copes(i).folder,'/');
    sub_id    = dir_split{7}; 
    copyfile(fullfile(Copes(i).folder,Copes(i).name), fullfile(Outdir, ['cope_' sub_id '.nii.gz'])); 
    command   = ['fslmaths ' fullfile(Varcopes(i).folder,Varcopes(i).name) ' -bin ' fullfile(Outdir, ['mask_' sub_id '.nii.gz'])];
    system(command);    
end

% Concatenate copes to 4D file
out_copes = cellstr(spm_select('FPList', Outdir, 'cope.*\.nii.gz$'));
command = strcat('fslmerge -t', {' '}, fullfile(Outdir,'copes.nii.gz'), {' '}, strjoin(out_copes'));
system(char(command));

% Create group mask
out_masks = cellstr(spm_select('FPList', Outdir, 'mask.*\.nii.gz$'));
command = strcat('fslmerge -t', {' '}, fullfile(Outdir,'group_mask.nii.gz'), {' '}, strjoin(out_masks'));
system(char(command));

command = ['fslmaths ' fullfile(Outdir,'group_mask.nii.gz') ' -Tmin ' fullfile(Outdir,'group_mask.nii.gz')];
system(command); 

HCP_contour_inf(Outdir,Outdir); 