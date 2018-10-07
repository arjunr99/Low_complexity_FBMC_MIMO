function [y] = vblast(x,H)

d = min(size(H));
y = zeros(d,1);
posn = zeros(1,d);

for l=1:d
    
    G = pinv(H);
    [~,pos] = min(diag(G*G'));
    posn(l) = pos;
    t = G(pos,:)*x;
    x = x - H(:,pos)*t;
    H = H(:,[1:pos-1 pos+1:end]);
    y(l) = t;
    
end

posn = vbreorder(posn);

y(posn) = y;
end