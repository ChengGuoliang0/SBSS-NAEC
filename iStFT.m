function Y = iStFT(X, N, win, shift)
nfft=(size(X,1)-1)*2;
L=size(X,2);
wlen=size(win,1);
W=zeros(N+wlen+nfft,1);
X=[X;conj(X(end-1:-1:2,:))];
Y=zeros(N+wlen+nfft,1);
for i=1:L
    sp=shift*i+1;
    tmp=real(ifft(X(:,i)));
    W(sp:sp+wlen-1)=W(sp:sp+wlen-1)+win;
    Y(sp:sp+wlen-1)=Y(sp:sp+wlen-1)+tmp(1:wlen);
end
W(shift*(L+1)+1:shift*(L+1)+wlen)=W(shift*(L+1)+1:shift*(L+1)+wlen)+win;
Y=Y(wlen+1:wlen+N)./W(wlen+1:wlen+N);
end