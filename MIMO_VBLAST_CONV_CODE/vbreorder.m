function y = vbreorder(x)

len = length(x);
pos = zeros(1,len);
y = zeros(1,len);

for l= 1:len
    
   ind = find(pos~=1);
   y(l) = ind(x(l));
   pos(y(l))= 1;    
   
end

end