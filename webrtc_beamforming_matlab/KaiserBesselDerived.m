function wint = KaiserBesselDerived(alpha, winlen)
%²úÉúKaiser Bessel Derived window

half = floor( ( winlen + 1 ) / 2 );
sum = 0.0;
wint = zeros(winlen,1);
for i = 1 : half + 1
    r = ( 4.0 * ( i - 1 ) ) / winlen - 1.0;
    sum = sum + I0( pi * alpha * real( sqrt( 1 - r ^ 2 ) ) );
    wint(i) = sum;
end

for i = winlen : -1 : half + 1
    wint( winlen - i + 1 ) = sqrt( wint( winlen - i + 1 ) / sum );
    wint( i ) = wint( winlen - i + 1 );
end

if( mod(winlen,2) == 1 )
    wint( half - 1 ) = sqrt( wint( half - 1 ) / sum );
end







