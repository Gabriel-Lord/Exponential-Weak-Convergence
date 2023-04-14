clear all
close all
%
T=1;
d=1;
m=d;
dx=1/(d+1);
dxsq=dx.*dx;
kappa=2;
DTMX=0.1;
ee=ones(d,1);
AA=0.5;
BB=0.1;
A = AA*(diag(ee(1:d-1),-1)-2*diag(ee)+diag(ee(1:d-1),1))./dxsq;
B = BB*eye(m); 
f = @(u) u-u.^3; %-u.^3;
g = @(u) u./(1+u.^2);

u0=rand(d,1);
x=dx:dx:1-dx;x=x';
u0=0.5*exp(-10*(x-0.5).^2);
%plot(u0);

epsilon=1.e-3;
%epsilon=1e-2;
%epsilon=5e-2;
%epsilon=7e-3;
%epsilon=5e-3;
				%
Levels=ceil(log(2*T/epsilon)/log(kappa))+1;
DT=T*kappa.^(-(0:Levels-1))
L0=find(DT<=DTMX, 1)

MINT=5000;

mrun=2;
for ii=2:-1:0
ii
 if ii==0
    EstimateWkErr=[];
    ETAM=[];
    EGBM=[];
    icase='MLMCL0'  % uses variance reduction on level 0 
  elseif  ii==1
    EstimateWkErr=[];
    ETAM=[];
    EGBM=[];
    icase='MLMC'    % MLMC for both tamed and GBM - done just once
  elseif  ii==2
    EstimateWkErr=[];
    ETAM=[];
    EGBM=[];
    icase='Trad.tamed' % MLMC for both tamed and GBM done for different epsilon values.
  end

  
  
  for kk=1:mrun

    if ii<2
      tic
      [EPuT, MT, S1T, S3T, Dt]=mlmc(u0,T,d,m,A,f,B,g,kappa,epsilon,DTMX,MINT,icase);
      trun(kk,1)=toc;
      True=EPuT(1);
    end
%
    if ii==2
%      MINT=1500;
      tic  % get true solution
      [EPuTB,MTB,S1TB,S3TB,DtB]=mlmc(u0,T,d,m,A,f,B,g,kappa,epsilon,DTMX,MINT,'Trad.tamed');
      trunTrue=toc;
      TrueB=EPuTB(1);
      for k=1:Levels-L0;
	EPS(k)=epsilon*pow2(k);
	[EPuTB,MTB,S1TB,S3TB,DtB]=mlmc(u0,T,d,m,A,f,B,g,kappa,EPS(k),DTMX,MINT,'Trad.gbm');
	DTB(k,1)=DtB(end);
%	L0=min(find(S3TB>0));
	estimate(k,kk)=EPuTB(2);
	EGBMB(k,kk)=abs(TrueB-estimate(k,kk));
      end
      trun(kk,1)=trunTrue+toc;
      tic
      for k=1:Levels-L0
	EPS(k)=epsilon*pow2(k);
	[EPuTB,MTB,S1TB,S3TB,DtB]=mlmc(u0,T,d,m,A,f,B,g,kappa,EPS(k),DTMX,MINT,'Trad.tamed');
	DTB(k)=DtB(end);
%	L0=min(find(S1TB>0));
	estimate(k,kk)=EPuTB(1);
	ETAMB(k,kk)=abs(TrueB-estimate(k,kk));
      end
      trun(kk,2)=trunTrue+toc;
%
%        
%
    elseif ii==1
      True=EPuT(1);
      for k=L0+2:length(Dt)
	EstimateWkErr(k,kk)=(True-(sum(S3T(L0:k)./MT(L0:k))));
	ETAM(k,kk)=sum(S1T(k:end)./MT(k:end));
	EGBM(k,kk)=sum(S3T(k:end)./MT(k:end));
      end
    elseif ii==0
      EstimateWkErr(:,kk)=True*ones(size(Dt));
      for k=L0+2:length(Dt)
	EstimateWkErr(k,kk)=(EstimateWkErr(k,kk)-(sum(S3T(1:k)./MT(1:k))));
      end
    end % if ii==
  end % for kk=1
  

  AvTime=mean(trun)
  StdTime=std(trun)
  

  if ii==2
    EEGBMB=mean(EGBMB,2);
    SGBMB=std(EGBMB');
    EETAMB=mean(ETAMB,2);
    STAMB=std(ETAMB');
    figure
    errorbar(DTB,EEGBMB,SGBMB,'b-','LineWidth',1.5)
    hold on
    errorbar(DTB,EETAMB,STAMB,'r-o','LineWidth',1.5)
    plot(DTB,DTB,'k')
    set(gca, 'XScale','log', 'YScale','log')
    hold off
    grid
    xlabel('Time step')
    ylabel('Weak Error')
    legend('GBM Tamed', 'Exponential Tamed','Reference slope 1')
    title([icase,': Av runtime=',num2str(AvTime),' (+/-',num2str(StdTime),')'])
    fname=['d',num2str(d),icase,'ii',num2str(ii),'MINT',num2str(MINT),'A',num2str(AA),'B',num2str(BB)];
    % save(['./',fname,'.mat'])
    saveas(gcf, ['./',fname,'.fig'], 'fig')
  else
    DT=Dt(L0+1:end);
    ERR=mean(abs(EstimateWkErr(L0+1:end,:)),2);
    STD=std(abs(EstimateWkErr(L0+1:end,:))');
    
    figure
    errorbar(DT,ERR,STD,'LineWidth',1.5)
    hold on
    plot(DT,DT,'k')
    set(gca, 'XScale','log', 'YScale','log')
    hold off
    grid
    xlabel('Time step')
    ylabel('Weak Error')
  
  legend(['Weak Error:',icase],'Reference slope 1')
  title([icase,'ii',num2str(ii),': Av runtime=',num2str(AvTime),' (+/-',num2str(StdTime),')'])
  fname=['d',num2str(d),icase,'ii',num2str(ii),'MINT',num2str(MINT),'A',num2str(AA),'B',num2str(BB)];
  saveas(gcf, ['./',fname,'.fig'], 'fig')

	end
  
  switch icase
    case 'MLMC'
      EETAM=mean(abs(ETAM),2);
      STDTAM=std(abs(ETAM)');
      EEGBM=mean(abs(EGBM),2);
      STDGBM=std(abs(EGBM)');

      figure
      errorbar(Dt,EEGBM(1:end),STDGBM(1:end),'b-','LineWidth',1.5)
      hold on
      errorbar(Dt,EETAM(1:end),STDTAM(1:end),'r-o','LineWidth',1.5)
      plot(DT,DT,'k')
      set(gca, 'XScale','log', 'YScale','log')
      hold off
      grid
      xlabel('Time step')
      ylabel('Weak Error')
      legend('GBM Tamed', 'Exponential Tamed','Reference slope 1')
%      title([icase,': Av runtime=',num2str(AvTime),' (+/-',num2str(StdTime),')'])
       fname=['d',num2str(d),icase,'SR','MINT',num2str(MINT),'A',num2str(AA),'B',num2str(BB)];
  saveas(gcf, ['./',fname,'.fig'], 'fig')
	     end
  
end

