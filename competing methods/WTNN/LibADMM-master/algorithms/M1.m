function  [error1,psnr1,error2,psnr2,itr1,itr]=M1(C1,C,omega)

 [X,obj,err,iter] =lrtc_tnn(C1,omega);
error1=norm(unfold(C-X),'fro')/norm(unfold(C),'fro');
psnr1=psnr(C,X);
itr1=iter;
 [X,obj,err,iter] =lrtcn_tnn(C1,omega,1);
error2(1)=norm(unfold(C-X),'fro')/norm(unfold(C),'fro');
psnr2(1)=psnr(C,X);
itr(1)=iter;
 [X,obj,err,iter] =lrtcn_tnn(C1,omega,2);
error2(2)=norm(unfold(C-X),'fro')/norm(unfold(C),'fro');
psnr2(2)=psnr(C,X);
itr(2)=iter;
 [X,obj,err,iter] =lrtcn_tnn(C1,omega,5);
error2(3)=norm(unfold(C-X),'fro')/norm(unfold(C),'fro');
psnr2(3)=psnr(C,X);
itr(3)=iter;
 [X,obj,err,iter] =lrtcn_tnn(C1,omega,10);
error2(4)=norm(unfold(C-X),'fro')/norm(unfold(C),'fro');
psnr2(4)=psnr(C,X);
itr(4)=iter;
 [X,obj,err,iter] =lrtcn_tnn(C1,omega,20);
error2(5)=norm(unfold(C-X),'fro')/norm(unfold(C),'fro');
psnr2(5)=psnr(C,X);
itr(5)=iter;

 [X,obj,err,iter] =lrtcn_tnn(C1,omega,50);
error2(6)=norm(unfold(C-X),'fro')/norm(unfold(C),'fro');
psnr2(6)=psnr(C,X);
itr(6)=iter;
disp(error1);
disp(psnr1);
disp(itr1);
disp(error2);
disp(psnr2);
disp(itr);


