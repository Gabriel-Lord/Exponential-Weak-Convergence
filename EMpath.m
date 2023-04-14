function [t, u, v]=EMpath(u0,T,N,d,m,A,fhandle,B,ghandle,kappa0,M,icase)
Dtref=T/N;   % small step
Dt=kappa0*Dtref;  % large step 
NN=N/kappa0;
u=zeros(d,M,NN+1); % Tamed
v=zeros(d,M,NN+1); % GBM
t=zeros(NN+1,1);
gdW=zeros(d,M);
sqrtDtref=sqrt(Dtref);
u_n=u0; v_n=u0;
Phi = expm((A-0.5*B^2)*Dt);
for n=1:NN+1
  t(n)=(n-1)*Dt;
  u(:,:,n)=u_n;
  v(:,:,n)=v_n;   
  dW=sqrtDtref*squeeze(sum(randn(m,M,kappa0),3));
  switch icase
    case {'Trad.tamed'}
      for mm=1:M
	gudW(:,mm)=(B*u_n(:,mm) + ghandle(u_n(:,mm)) ).*dW(:,mm);    
      end
      fu=fhandle(u_n);
      tamedfu=fu./(1+Dt*norm(fu));
      u_new=u_n+Dt*(A*u_n + tamedfu) + gudW;
      u_n=u_new;
    case {'Trad.gbm'}
      fv=fhandle(v_n);
      tamedfv=fv./(1+Dt*norm(fv));
      for mm=1:M
	gv=ghandle(v_n(:,mm));
	gvdW=gv.*dW(:,mm);
	v_new(:,mm)=Phi*expm(diag(B*dW(:,mm)))*(v_n(:,mm)+Dt*(tamedfv(:,mm)-B*gv)+gvdW);
      end
      v_n=v_new;
    otherwise
      fv=fhandle(v_n);
      tamedfv=fv./(1+Dt*norm(fv));      
      for mm=1:M
	gudW(:,mm)=(B*u_n(:,mm) + ghandle(u_n(:,mm)) ).*dW(:,mm);    
	gv=ghandle(v_n(:,mm));
	gvdW=gv.*dW(:,mm);
	v_new(:,mm)=Phi*expm(diag(B*dW(:,mm)))*(v_n(:,mm)+Dt*(tamedfv(:,mm)-B*gv)+gvdW);
      end
      fu=fhandle(u_n);
      tamedfu=fu./(1+Dt*norm(fu));
      u_new=u_n+Dt*(A*u_n + tamedfu) + gudW;
      u_n=u_new;
      v_n=v_new;
  end
end

