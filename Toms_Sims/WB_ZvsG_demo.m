% Demo of difference between performance of Gaussian (Z) and Rademacher 
% Wild (Multiplier) Bootstrap
%
% Setting:
% nSub subjects each with a nVox-voxel image with independent voxels,
% with nB bootstrap samples.  Entire enterprise is repeated nMC times.
%
% Goal: Estimate upper quantiles of the maximum distribution of one-sample t 
% statistic
%
% Results: Use both true, known maximum distribution, Monte Carlo and
% Bootstrap samples to estimate 0.8, 0.9 and 0.95 quantiles.
%

clear
close all
set(0, 'DefaultFigureRenderer', 'painters');


nSub=60;
nVox=1000;
nB=1000;
nMC=500;
Alph=[0.8 0.9 0.95];

Tstat = @(X)mean(X)./std(X)*sqrt(size(X,1));
Xbar  = @(X)mean(X)*sqrt(size(X,1));

Tmax=@(X,df,nV)tcdf(X,df).^nV;

Y=randn(nSub,nVox,nMC);
[MTG,MTR]=deal(zeros(nB,nMC));


for i=1:nB

  % Standardised residuals
  res = (Y-mean(Y))./std(Y);

  if i==1

    MTG(1,:)=max(Tstat(Y));
    MTR(1,:)=MTG(1,:);

  else

    YG=Y.*   randn([nSub 1 nMC]);
    YR=Y.*(2*randi(2,[nSub 1 nMC])-3);

    MTG(i,:)=max(Tstat(YG));
    MTR(i,:)=max(Tstat(YR));

  end

  fprintf('.')

end

% Find true critical values

MnMx=[min(MTG(:)), max(MTG(:))];

optopt = optimoptions('fmincon','Display','off');
u80 = fmincon(@(x)(Tmax(x,nSub-1,nVox)-0.80).^2,mean(MnMx),[],[],[],[],MnMx(1),MnMx(2),[],optopt);
u90 = fmincon(@(x)(Tmax(x,nSub-1,nVox)-0.90).^2,mean(MnMx),[],[],[],[],MnMx(1),MnMx(2),[],optopt);
u95 = fmincon(@(x)(Tmax(x,nSub-1,nVox)-0.95).^2,mean(MnMx),[],[],[],[],MnMx(1),MnMx(2),[],optopt);


xx=linspace(MnMx(1),MnMx(2),100);
xMC=(1:nMC)/(nMC+1);
xB=(1:nB)/(nB+1);


set(gcf,'PaperSize',[4 3]*2.54,'PaperPosition',[0 0 [4 3]*2.54])

h1=plot(sort(MTG),xB);title('Gaussian')
hold on;
plot(prctile(MTG,95),0.95,'.','markersize',10)
h=[abline('h',[0.8 0.9 0.95]);abline('v',[u80 u90 u95])];
%set(h,'linestyle','-');
uistack(h,'top')
h2=plot(xx,Tmax(xx,nSub-1,nVox),'color',[0.8 0 0 ],'linewidth',4);
h3=plot(sort(MTG(1,:)),xMC,'--','color',[0.5 1 0.5],'linewidth',3);
ylim([0.75 1]);xlim(MnMx)
legend([h2 h3 h1(1)],'F_{Max} true','F_{Max} MC','F_{Max} Boot','location','southeast')
hold off
print -dpdf WB_Gauss.pdf

h1=plot(sort(MTR),xB);title('Rademacher')
hold on;
plot(prctile(MTR,95),0.95,'.','markersize',10)
h=[abline('h',[0.8 0.9 0.95]);abline('v',[u80 u90 u95])];
%set(h,'linestyle','-');
uistack(h,'top')
h2=plot(xx,Tmax(xx,nSub-1,nVox),'color',[0.8 0 0 ],'linewidth',4);
h3=plot(sort(MTG(1,:)),xMC,'--','color',[0.5 1 0.5],'linewidth',3);
ylim([0.75 1]);xlim(MnMx)
legend([h2 h3 h1(1)],'F_{Max} true','F_{Max} MC','F_{Max} Boot','location','southeast')
hold off

print -dpdf WB_Rade.pdf
