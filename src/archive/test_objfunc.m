function obj=test_objfunc(x,timevec,resvec,simui)
  simusig=zeros(size(timevec));
  for i=1:simui
    shifti=(i-1)*3;
    simusig=simusig+(x(3+shifti).*sin(2*pi*x(1+shifti)*timevec)).*exp(-x(2+shifti)*timevec);
  end
  obj= -log(sum((simusig-resvec).^2));
