function [e, W] = SBSS_NAEC(x, y, p, nfft, eta)

% Semi-blind source separation algorithm for nonlinear acoustic echo cancellation
%
%   e : estimate of the near-end signal (L x 1)
%   W : demixing matrices (p+1 x p+1 x nfft/2+1)
%   x : far-end input signal (L x 1)
%   y : microphone signal (L x 1)
%   p : expansion order
%   nfft : number of FFT points
%   eta : learning rate

L=length(x);
xp=zeros(L,p);
for i=1:p
    xp(:,i)=x.^(2*i-1);
end
yxp=[y,xp];

% Short-time Fourier transform
win=hanning(nfft,'periodic');
shift=fix(nfft/4);
N=fix((L+nfft)/shift)-1;
nf=fix(nfft/2)+1;
Y=zeros(p+1,N,nf);
for i=1:p+1
    Y(i,:,:)=StFT(yxp(:,i),nfft,win,shift).';
end

W=zeros(p+1,p+1,nf);
En=zeros(p+1,1,nf);
R=zeros(size(W));
dW=zeros(size(W));
E=zeros(nf,N);

% Initialize
for k=1:nf
    W(:,:,k)=eye(p+1);
end

% Online SBSS algorithm
for n=1:N
    Yn=Y(:,n,:);
    for k=1:nf
        En(:,:,k)=W(:,:,k)*Yn(:,:,k);
    end
    Ssq=sqrt(sum(abs(En).^2,3));
    Ssq1=(Ssq+1e-6).^-1;
    for k=1:nf
        
        % Compute multivariate score function
        Phi=Ssq1.*En(:,:,k);
        
        % Compute scaling factors
        R(:,:,k)=Phi*En(:,:,k)';
        dk=sum(sum(abs(R(:,:,k))))/(p+1);
        ck=1/dk;
        
        % Update demixing matrices using constrained scaled natural gradient strategy
        dW(:,:,k)=(eye(p+1)-R(:,:,k)/dk)*W(:,:,k);
        dW(2:p+1,:,k)=0;
        W(:,:,k)=ck*(W(:,:,k)+eta*dW(:,:,k));
        W(1,:,k)=W(1,:,k)/W(1,1,k);
        W(2:p+1,2:p+1,k)=eye(p);
        
        E(k,n)=W(1,:,k)*Yn(:,:,k);
    end
end

% Reconstruct the estimated near-end signal using inverse STFT
e=iStFT(E,L,win,shift);

end