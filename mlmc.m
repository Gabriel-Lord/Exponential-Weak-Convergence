function [EPuv, M, S1, S3, DT]=mlmc(u0,T,d,m,A,fhandle,B,ghandle,kappa,epsilon,DTMX,MINT,icase)
%
  Levels=ceil(log(2*T/epsilon)/log(kappa))+1;
  DT=T*kappa.^(-(0:Levels-1)); DT=DT';
  CT=sqrt(T*(1+1./kappa)); CTI=1./CT;
  L0=find(DT<=DTMX, 1);
  S1=zeros(Levels,1); S2=zeros(Levels,1);
  S3=zeros(Levels,1); S4=zeros(Levels,1);

  M=MINT*ones(Levels,1); ML=zeros(Levels,1); VL=zeros(Levels,1);
				%
  for L=L0:Levels
    N=kappa^L;   
    [S1,S2,S3,S4]=getlevel(u0,T,N,d,m,A,fhandle,B,ghandle,kappa,M(L),L,L0,S1,S2,S3,S4,icase); 
    TMP12=S2(L0:L)./M(L0:L)-(S1(L0:L)./M(L0:L)).^2;
    TMP34=S4(L0:L)./M(L0:L)-(S3(L0:L)./M(L0:L)).^2;
    VL(L0:L)=max(TMP12,TMP34);
    KC=2*(sqrt(VL(L0)*T/DT(L0)))/epsilon^2;

    if L>L0
      KC=KC+(2*CTI*sum(sqrt(VL(L0+1:L)./DT(L0+1:L))))/epsilon^2;
      ML(L0+1:L)=ceil(KC*CT*sqrt(VL(L0+1:L).*DT(L0+1:L)));
      %ML(L0+1:L)=ceil(5*KC*CT*sqrt(VL(L0+1:L).*DT(L0+1:L)));
    end				%
    ML(L0)=ceil(KC*sqrt(VL(L0)*DT(L0)/T));

    for l=L0:L
      dM=ML(l)-M(l);
      if dM>0
	N=kappa^l;   M(l)=M(l)+dM;
	[S1,S2,S3,S4]=getlevel(u0,T,N,d,m,A,fhandle,B,ghandle,kappa,dM,l,L0,S1,S2,S3,S4,icase);
      end
    end
    EPuv(1)=sum(S1(L0:L)./M(L0:L));
    EPuv(2)=sum(S3(L0:L)./M(L0:L));
%    [L EPuv]
  end
