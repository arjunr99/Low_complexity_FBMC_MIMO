function [y] = vblast_rec(x,H)

d = min(size(H));

if(d>1)
    
    G = pinv(H);
    [~,pos] = min(diag(G*G'));
    t = G(pos,:)*x;
    x1 = x - H(:,pos)*t;
    H1 = H(:,[1:pos-1 pos+1:end]);
    temp = vblast_rec(x1,H1);
    y = [t temp];
else
    y = H'*x/norm(H)^2;
end

end