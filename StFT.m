function Y = StFT(X, nfft, win, shift)
wlen=size(win,1);
N=size(X,1);
X=[zeros(wlen,1);X;zeros(nfft,1)];
L=fix((N+wlen)/shift)-1;
Y=zeros(nfft,L);
Xi=zeros(nfft,1);
for i=1:L
    sp=shift*i+1;
    Xi(1:wlen)=win.*X(sp:sp+wlen-1);
    Y(:,i)=fft(Xi);
end
Y=Y(1:fix(nfft/2)+1,:);
end