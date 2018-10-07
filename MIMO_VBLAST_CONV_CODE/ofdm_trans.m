function ofdm_sig=ofdm_trans(input , fft_size , cpLength)

inSize=length(input);

if(mod(inSize,fft_size))
    
    input=[input zeros(1,fft_size-mod(inSize,fft_size))];

end

inSize=length(input);

NumOfdmSym=inSize / fft_size;
LenOfdmSym = fft_size + cpLength;

ofdm_sig=zeros(1,inSize + (NumOfdmSym * cpLength));

for l=1:NumOfdmSym
    
    ofdm_sig(LenOfdmSym * (l-1) + cpLength + 1 : LenOfdmSym * l)=ifft(...
        input(fft_size * (l-1) + 1 : fft_size * l));
    ofdm_sig(LenOfdmSym * (l-1) + 1 : LenOfdmSym * (l-1) + cpLength) = ...
        ofdm_sig(LenOfdmSym * l - cpLength + 1 : LenOfdmSym * l);
    
end

end