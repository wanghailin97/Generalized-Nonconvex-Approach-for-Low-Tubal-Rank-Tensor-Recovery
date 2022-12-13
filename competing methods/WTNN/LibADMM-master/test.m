clear
A=rand(1000,200,20);
B=rand(200,100,20);
C=tprodg(A,B);
tic;[X, tnn, trank] = prox_tnn(C,0.05),timee1=toc;
tic;[X2, tnn2, trank2] = prox_tnng(C,20),timee2=toc;