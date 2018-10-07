function F = init(k)

F(1) = sum(k(1:31))  + 0.5;
for i = 1:16
F(i+1) = k(i)^2 + k(32-i)^2 - 1;
end