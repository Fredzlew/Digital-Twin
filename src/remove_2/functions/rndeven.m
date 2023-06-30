function y = rndeven(x)
    x = floor(x);
    x(x <= 1) = 2;
    y = mod(x,2)+x;
end