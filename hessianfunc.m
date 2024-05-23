function yprpr = hessianfunc(x)
    syms x1 x2
    F = 100 * (x2 - x1^2)^2 + (6.4 * (x2 - 0.5)^2 - x1 - 0.6)^2;
    H = hessian(F);
    yprpr = double(subs(H, [x1 x2], [x(1) x(2)]));
end
