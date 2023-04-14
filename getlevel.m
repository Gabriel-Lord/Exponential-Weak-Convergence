function [S1,S2,S3,S4]=getlevel(u00,T,N,d,m,A,fhandle,B,ghandle,...
                              kappa,MS,L,L0,S1,S2,S3,S4,icase)

  SA=0; SB=0; SC=0; SD=0;
  Mstep=10000;
for ktmp=1:floor(MS/Mstep)+1;
%  parfor ktmp=1:floor(MS/Mstep)+1;
    M=1+(ktmp-1)*Mstep;
    MM=min(Mstep,MS-M+1);
    u0=bsxfun(@times,u00,ones(d,MM)); % initial data
    if L==L0
      [t,u,v]=EMpath(u0, T, N, d, m,A,fhandle, B ,ghandle,1,MM,icase);
      u=squeeze(u(:,:,end));
      v=squeeze(v(:,:,end));
      switch icase
	case 'MLMCL0'
	  SA=SA+sum(phi(u)-phi(v));
	  SB=SB+sum((phi(u)-phi(v)).^2);
	case 'MLMC'
	  SA=SA+sum(phi(u));
	  SB=SB+sum((phi(u)).^2);
	  SC=SC+sum(phi(v));
	  SD=SD+sum(phi(v).^2);
	case 'Trad.tamed'
	  SA=SA+sum(phi(u));
	  SB=SB+sum((phi(u)).^2);
	case 'Trad.gbm'
	  SC=SC+sum(phi(v));
	  SD=SD+sum(phi(v).^2);
      end
    else
      defaultStream = RandStream.getGlobalStream;
      savedState = defaultStream.State;
      [t,uu,vv]=EMpath(u0, T, N, d, m,A, fhandle,B,ghandle,1,MM,icase);
      uref=squeeze(uu(:,:,end));
      vref=squeeze(vv(:,:,end));
      defaultStream.State = savedState;
      [t,uu,vv]=EMpath(u0, T, N, d, m,A,fhandle,B,ghandle,kappa,MM,icase);
      u=squeeze(uu(:,:,end));
      v=squeeze(vv(:,:,end));

      Xu=(phi(uref)-phi(u));
      Xv=(phi(vref)-phi(v));
      SA= SA+sum(Xu);
      SB= SB+sum(Xu.^2);
      SC= SC+sum(Xv);
      SD= SD+sum(Xv.^2);
    end
  end  
  S1(L)=S1(L)+SA;
  S2(L)=S2(L)+SB;
  S3(L)=S3(L)+SC;
  S4(L)=S4(L)+SD;  
  
  function phiv=phi(v)
    phiv=sqrt(diag(v'*v));
