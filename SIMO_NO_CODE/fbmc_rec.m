%%%% change shift_pts value to change spectral efficiency


function [rec] = fbmc_rec(received,Ncarriers,fb)

%len=length(received);

%lf = size(fb,1);

%shift_pts = 0 for 32/32, 1 for 32/33, 2 for 32/34
shift_pts = 1;

up = (32+shift_pts)*Ncarriers/32;

fb = fliplr(fb');

for l=1:Ncarriers
         
    temp = conv(received, fb(l,:),'valid');
    b1 = downsample( temp , up);
    rec_sig( : , l )= b1.';
    
end

rec = rec_sig.';

rec=rec(:).';
end
