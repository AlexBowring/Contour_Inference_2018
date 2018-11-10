clear
close all
set(0, 'DefaultFigureRenderer', 'painters');
maxNumCompThreads(8);

nSub = 60;
nVox = 1000;
nB   = 1000;
nMC  = 1000;

% This is where we started
Resid=true;
JASA=true;
WB_ZvsG(nSub,nVox,nB,nMC,Resid,JASA)

% This is what Fabian proposed
Resid=true;
JASA=false;
WB_ZvsG(nSub,nVox,nB,nMC,Resid,JASA)

% This makes WB work optimally
Resid=false;
JASA=false;
WB_ZvsG(nSub,nVox,nB,nMC,Resid,JASA)

