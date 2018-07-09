function a = gcd( a, b )
%UNTITLED Summary of this function goes here
while ( b )
    temp = a;
    a = b;
    b = mod( temp, b );
end

end

