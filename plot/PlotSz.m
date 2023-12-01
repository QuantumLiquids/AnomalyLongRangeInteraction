L = 128;
sz = zeros(1, L);
for i = 1:L
 sz(i) = data{i}{2}(1);
end
plot(sz,'-o');