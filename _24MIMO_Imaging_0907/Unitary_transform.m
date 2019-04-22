function Z = Unitary_transform(Y)
[P,~] = size(Y);
if mod(P,2) == 0
    m = P/2;
else
    m = (P-1)/2;
    y = Y(m+1,:);
end
Y1 = Y(1:m,:);
Y2 = Y((end-m+1):end,:);
Y21 = IEM(m) * Y2;
if mod(P,2) == 0
    Z = [real(Y1+Y21) -imag(Y1+Y21);
        imag(Y1-Y21) real(Y1-Y21)];
else
    Z = [real(Y1+Y21) -imag(Y1+Y21);
        sqrt(2)*real(y) -sqrt(2)*image(y);
        imag(Y1-Y21) real(Y1-Y21)];
end