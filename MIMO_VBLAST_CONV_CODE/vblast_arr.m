function [y] = vblast_arr(x,H)

s = size(x);

y = zeros(min(size(H)),s(2));

for l=1:s(2)
    
    y(:,l) = vblast(x(:,l),H);

end

