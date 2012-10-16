function y = symlog(x)

neg = x < 0;
z = x == 0;
y = log10(abs(x)+1);
y(z) = 0;
y(neg) = -y(neg);