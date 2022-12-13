function OUT = GenLowTuba(n1,n2,n3,r)
P = normrnd(0,1/sqrt(n1),int32([n1,r,n3]));
Q = normrnd(0,1/sqrt(n1),int32([r,n2,n3]));
P = fft(P,[],3);
Q = fft(Q,[],3);
C = zeros(n1,n2,n3);
for i = 1 : n3
    C(:,:,i) = P(:,:,i)*Q(:,:,i);
end
OUT = ifft(C,[],3);


