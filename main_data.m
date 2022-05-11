num = 5;
imps = zeros(1, num);
for i=1:num
    imps(i) = realistic_frequencies();
end
imps

% mult 1 [2.2336    1.1489    1.5633    5.0113    3.4929] tau 5e-4
% mult 10 [0.9777    1.8556    1.0693    0.5650    0.2224] tau 20e-4