function ypr = gradfunc(x)
    syms x1 x2
    F = 100 * (x2 - x1^2)^2 + (6.4 * (x2 - 0.5)^2 - x1 - 0.6)^2;
    g = gradient(F);
    ypr = double(subs(g, [x1 x2], [x(1) x(2)]));
end
