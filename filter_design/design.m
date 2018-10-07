%function k = design
close all;
format long;
%TF=33/32 design
% k1-k31 anlysis 
Ncarriers = 128;
over = 32;
N = Ncarriers*33/32;
x0 = [-0.99 ,0.98, -0.97 ,0.96, -0.95 ,0.94, -0.93, 0.92,-0.91,0.9,-0.89,0.86,-0.82,0.78,-0.74,0.7,-0.65,0.6,-0.55,0.5,-0.45,0.4,-0.35,0.3,-0.25,0.2,-0.16,0.12,-0.09,0.06,-0.03];
%x0 = [-1 ,1, -1 ,1, -1 ,1, -1, 1,-1,1,-1,1,-1,1,-1,0.7,-07,0.01,-0.01,0.01,-0.01,0.01,-0.01,0.01,-0.01,0.01,-0.01,0.01,-0.01,0.01,-0.01];
%x0=randn(1,32);
options = optimoptions('fsolve','Display','iter');
options.MaxFunEvals = 20060;  % iteration
options.MaxIter = 50000;      % iteration
options.TolFun = 1e-12;       % tolerance
x=fsolve(@init,x0,options);
k0 = [x];
fun = @solv;
options.TolX = 1e-12;

k = fsolve(fun,k0,options);

%k = [-0.995200546670935,0.980826733816463,-0.957039977049698,0.924049099777236,-0.882196084881446,0.831866870263224,-0.773572288358335,0.707863959994013,-0.635404055460667,0.556903206164610,-0.473159181437296,0.385073519500568,-0.293626859094314,0.200344289692203,-0.108190662743724,0.0329053673901957,-0.00128513833997250,-0.000160873211954348,-2.40257471908938e-05,-8.04692024328481e-06,4.00178835746455e-06,1.23156900112914e-06,5.43207467362760e-06,2.14042551829416e-06,4.84979918876288e-06,2.02092920188864e-06,4.23144863933846e-06,1.79664220352496e-06,3.79683971241770e-06,1.63229474467445e-06,3.55950013353687e-06];

H1 =[1 k(1:31) zeros(1,(N-2)*over+1) k(31:-1:1)];
%H2 =[1 k(32:62) zeros(1,(N-2)*over+1) k(62:-1:32)];
y1 = ifft(H1);
%y2 = ifft(H2);

f_l = 10;
fft_size = f_l*N*over;

Y1 = fft(y1,fft_size);
%Y2 = fft(y2,fft_size);

cumul =zeros(size(Y1));
for l=1:N
   cumul = cumul + circshift(Y1.*conj(Y1),[0,l*over*f_l]); 
end

%plot(abs(cumul));

%plot(corr);
corr = conv(y1,fliplr(y1));
y=downsample(circshift(corr,[0,-(over*N-1)]),N);
%F(4) = sum(abs(y(2:end)));

H3 = circshift(H1,[0,over+1]);
y3 = ifft(H3);
corr2 = conv(y1,conj(fliplr(y3)));
y4=downsample(circshift(corr2,[0,-(over*N-1)]),N);

2*sum(abs((y(2:over))/(y1*y1')))   % ISI
sum(abs((y4(1:over)/(y1*y1'))))    % ICI
%k;

%  save('filter_final.mat');
%  exit;