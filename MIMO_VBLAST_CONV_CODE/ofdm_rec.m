function y=ofdm_rec(input,fft_size,cpLength)

symLen=fft_size + cpLength;

numSym=length(input) / symLen; 

y=zeros(1 , fft_size * numSym);

for l=1:numSym
    
    y(fft_size * (l-1) + 1:fft_size * l)=fft(...
        input(symLen * (l-1) + cpLength + 1 : symLen * l));
    
end

end