%%%% to change the spectral efficiencies, change the value of shift_pts


function [sent_sig,s]=fbmc_trans( x , Ncarriers,fb)

sent = x.';

%shift_pts = 0 for 32/32, 1 for 32/33, 2 for 32/34
shift_pts = 1;

up = (32+shift_pts)*Ncarriers/32;
se = size(sent);

if(se(1)>1)
    s = upsample(sent,up);
    
else
    s = sent;
end

sent_sig = zeros(length(s(:,1)) + length(fb(:,1))-1,1);

for l=1:Ncarriers
     sent_sig = sent_sig + conv(s(:,l),fb(:,l));
end

sent_sig=sent_sig.';

end
