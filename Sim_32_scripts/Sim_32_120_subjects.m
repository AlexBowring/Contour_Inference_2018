nSubj=120;
nRlz=30;

tic
cd /storage/maullz/Contour_Inference_2018/
addpath('/storage/essicd/spm8/')
addpath(genpath('/storage/maullz/Contour_Inference_2018/'))

rng('shuffle')
tID=str2num(getenv('SGE_TASK_ID'))
Sim_32(nSubj,['Sim_32_' num2str(nSubj) '_subjects_' sprintf('%03d',tID)],nRlz)
toc
