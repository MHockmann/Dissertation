function B = padarray(A,padsize)
%padarry Pads an array by zeros 
%
[m,n]=size(A);
A=[zeros(m,padsize) A zeros(m,padsize)];
B=[zeros(padsize,n+2*padsize); A;zeros(padsize,n+2*padsize)];
end

