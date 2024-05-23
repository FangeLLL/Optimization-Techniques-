function y = func(x)
    x1 = x(1);
    x2 = x(2);
    y = 100 * (x2 - x1^2)^2 + (6.4 * (x2 - 0.5)^2 - x1 - 0.6)^2;
end
