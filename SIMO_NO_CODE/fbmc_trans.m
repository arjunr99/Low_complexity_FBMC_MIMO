%%%% change shift_pts value to change spectral efficiency


function [sent_sig,s]=fbmc_trans( x , Ncarriers,fb)

len_x = length(x);

if(mod(len_x,Ncarriers))
    x=[x zeros(1,Ncarriers-mod(len_x,Ncarriers))];    
end

len_x=length(x);

sent = (reshape(x,Ncarriers,len_x/Ncarriers)).';


%lf = size(fb,1);
%shift_pts = 0 for 32/32, 1 for 32/33, 2 for 32/34
shift_pts = 1;

up = (32+shift_pts)*Ncarriers/32;
  
if((len_x/Ncarriers) > 1)
%sent_c = sent;
    s=upsample(sent,up);
else
    s = sent;
end

s_size = size(s);
sent_sig = zeros(length(s(:,1)) + length(fb(:,1))-1,1);

for l=1:Ncarriers
     sent_sig = sent_sig + conv(s(:,l),fb(:,l));
end

sent_sig=sent_sig.';

end
