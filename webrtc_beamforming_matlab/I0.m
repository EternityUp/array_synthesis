function y = I0(x)
%Modified Bessel function of order 0 for complex inputs

y = x / 3.75;
y = y * y;
y = 1.0 + y * ( ...
    3.5156229 + y * ( ...
    3.0899424 + y * ( ...
    1.2067492 + y * ( ... 
    0.2659732 + y * ( ...
    0.360768e-1 + y * 0.45813e-2)))));



