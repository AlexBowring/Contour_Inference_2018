clear
close all
set(0, 'DefaultFigureRenderer', 'painters');
maxNumCompThreads(8);

nSub = 60;
nVox = 1000;
nB   = 1000;
nMC  = 50;

% This is where we started
Resid=true;
JASA=true;
WB_ZvsG_demo

% This is what Fabian proposed
Resid=true;
JASA=false;
WB_ZvsG_demo

% This makes WB work optimally
Resid=false;
JASA=false;
WB_ZvsG_demo

