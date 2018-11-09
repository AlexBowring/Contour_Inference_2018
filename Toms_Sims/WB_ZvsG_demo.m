function WB_ZvsG(nSub,nVox,nB,nMC,Resid,JASA)
% Demo of difference between performance of Gaussian (Z) and Rademacher 
% Wild (Multiplier) Bootstrap
%
% Inputs
%   nSub  - Number of subjects (e.g. 60)
%   nVox  - Number of 'voxels' (independent tests) to search over (e.g. 1000)
%   nB    - Number of Bootstrap samples (e.g. 9999)
%   Resid - Should raw data be residualised
%   JASA  - Should the statistic be as in JASA paper, or as per Fabian's proposal
%
% Setting:
%   nSub subjects each with a nVox-voxel image with independent voxels,
%   with nB bootstrap samples.  Entire enterprise is repeated nMC times.
%
% Goal:
%   Estimate upper quantiles of the maximum distribution of one-sample t statistic
%
% Results:
%   Use both true, known maximum distribution, Monte Carlo and Bootstrap samples
%   to estimate 0.8, 0.9 and 0.95 quantiles. 
%
% If JASA==true the WB samples are summed and scaled by n^(-1/2) (assuming true variance 1.0)
% If JASA==false the WB samples are averaged, divided by stdev and scaled by n^(1/2) (reflecting 
% actual studentisation in confidence statement)

Alph = [0.8 0.9 0.95];

if Resid, str1='Resid'; else str1='Raw'; end
if JASA;  str2='JASA';  else str2='Stud'; end
fnm=sprintf('S%d_v%d_B%d_MC%d_%s_%s',nSub,nVox,nB,nMC,str1,str2);

if JASA
  Stat = @(X)sum(X)/sqrt(size(X,1));
else
  Stat = @(X)mean(X)./std(X)*sqrt(size(X,1));
end

MaxCDF=@(X,df,nV)tcdf(X,df).^nV;

% Monster data matrix
% Subjects x Voxels x Monte Carlo realisations
Y=randn(nSub,nVox,nMC);

% Outputs: Max stat MC & Bootstrap G & R
MTMC = zeros(1,nMC);
[MTG,MTR]=deal(zeros(nB,nMC));

if Resid
  % If we center the original data for Monte Carlo, the statistic will be 
  % exactly zero, so we just scale
  MTMC = deal(max(Stat(Y./std(Y))));
  % Now standardise for good
  Y = (Y-mean(Y))./std(Y);
else
  MTMC = deal(max(Stat(Y)));
end

for i=1:nB

  YG = Y .*    randn(  [nSub 1 nMC]);
  YR = Y .* (2*randi(2,[nSub 1 nMC])-3);

  MTG(i,:) = max(Stat(YG));
  MTR(i,:) = max(Stat(YR));
  
  fprintf('.')

end

%
% Plotting
%

MnMx=[min(MTG(:)), max(MTG(:))];

% Find true critical values
if JASA
  df = Inf;
else
  df = nSub-1;
end

optopt = optimoptions('fmincon','Display','off');
u80 = fmincon(@(x)(MaxCDF(x,df,nVox)-0.80).^2,mean(MnMx),[],[],[],[],MnMx(1),MnMx(2),[],optopt);
u90 = fmincon(@(x)(MaxCDF(x,df,nVox)-0.90).^2,mean(MnMx),[],[],[],[],MnMx(1),MnMx(2),[],optopt);
u95 = fmincon(@(x)(MaxCDF(x,df,nVox)-0.95).^2,mean(MnMx),[],[],[],[],MnMx(1),MnMx(2),[],optopt);

xx   = linspace(MnMx(1),MnMx(2),100);
xxMC = (1:nMC)/(nMC+1);
xxB  = (1:nB)/(nB+1);

set(gcf,'PaperSize',[4 3]*2.54,'PaperPosition',[0 0 [4 3]*2.54])

h1=plot(sort(MTG),xxB);title('Gaussian')
hold on;
plot(prctile(MTG,95),0.95,'.','markersize',10)
h=[abline('h',[0.8 0.9 0.95]);abline('v',[u80 u90 u95])];
%set(h,'linestyle','-');
uistack(h,'top')
h2=plot(xx,MaxCDF(xx,df,nVox),'color',[0.8 0 0 ],'linewidth',2);
h3=plot(sort(MTMC),xxMC,'--','color',[0.5 1 0.5],'linewidth',2);
ylim([0.75 1]);xlim(MnMx)
legend([h2 h3 h1(1)],'F_{Max} true','F_{Max} MC','F_{Max} Boot','location','southeast')
hold off
print('-dpdf',['WB_Gauss_' fnm '.pdf']);

h1=plot(sort(MTR),xxB);title('Rademacher')
hold on;
plot(prctile(MTR,95),0.95,'.','markersize',10)
h=[abline('h',[0.8 0.9 0.95]);abline('v',[u80 u90 u95])];
%set(h,'linestyle','-');
uistack(h,'top')
h2=plot(xx,MaxCDF(xx,df,nVox),'color',[0.8 0 0 ],'linewidth',2);
h3=plot(sort(MTMC),xxMC,'--','color',[0.5 1 0.5],'linewidth',2);
ylim([0.75 1]);xlim(MnMx)
legend([h2 h3 h1(1)],'F_{Max} true','F_{Max} MC','F_{Max} Boot','location','southeast')
hold off

print('-dpdf',['WB_Rade_' fnm '.pdf']);
