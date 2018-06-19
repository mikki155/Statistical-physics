%Ising model on ring
T = 0.2:0.032:4;
n = 2;%Must be between 2 and 5 because of the way the code is implemented!
h = 0;
J = 1;
H = zeros(1, 2^(n^2));
Z = zeros(1, length(T));
He = zeros(1, length(T));
C = zeros(1, length(T));

M = (dec2bin(0:(2^(n^2))-1)=='1');%Solves the configuration of -1 and 1 problem
M = double(M);
M(M == 0) = -1;

for i = 1:(2^(n^2))
    for j = 1:n^2
        if j < n^2
            H(i) = H(i) + M(i, j)*M(i, j+1);
        else
            H(i) = H(i) + M(i, n^2)*M(i, 1);
        end
    end
    Z = Z + exp(-H(i)./T);
    He = He + H(i).*exp(-H(i)./T);
    C = C + (H(i)^2.*exp(-H(i)./T));
end

C = (C./Z - (He./Z).^2)./(T.^2);

plot(T, C)
title('Heat capacity');
xlabel('T');
ylabel('C');
hold on