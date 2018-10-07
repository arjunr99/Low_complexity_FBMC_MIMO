function F= solv(k)
Ncarriers = 128;
over = 32;
N = Ncarriers*33/32;
for i = 1:16
F(i) = k(i)*k(i+31)+ k(32-i)*k(63-i)-1;
end
F(1)=sum(k(1:31))+0.5;
F(2)=sum(k(32:62))+0.5;

H1 =[1 k(1:31) zeros(1,(N-2)*over+1) k(31:-1:1)];

y1 = ifft(H1);

corr = conv(y1,fliplr(y1));
y=downsample(circshift(corr,[0,-(over*N-1)]),N);   % autocorrelation downsampled version

H3 = circshift(H1,[0,over+1]);
y3 = ifft(H3);
corr2 = conv(y1,conj(fliplr(y3)));
y4=downsample(circshift(corr2,[0,-(over*N-1)]),N);  % crosscorrealtion
%alpha =0.5;

F(1) = norm(y(3:over)) ; %+norm(y4(1:over)); %ISI
F(2) = norm(y4(1:over)); %ICI

end
