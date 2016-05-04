function cMap = interpMap(colorStart, colorEnd, n)

for i = 1:3
    cMap(1:n,i) = linspace(colorStart(i), colorEnd(i), n);
end
end