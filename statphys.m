%Partitioning function for two-dimensional lattice
T = 0.1:0.032:8;
n = 3;
h = 0.5;
J = 1;
if n > 2
    J = J/2;%This is because I calculated the spins twice for n>2
end
H = zeros(1, 2^(n^2));
Z = zeros(1, length(T));
He = zeros(1, length(T));
C = zeros(1, length(T));
m = zeros(1, length(T));

M = (dec2bin(0:(2^(n^2))-1)=='1');%Solves the configuration of -1 and 1 problem
M = double(M);
M(M == 0) = -1;

for i = 1:(2^(n^2))
    for j = 1:(n^2)
       if n == 2
          if j == 1
             H(i) = H(i) - J*M(i, j)*(M(i, j+1) + M(i, j+2));
          elseif j == 2
             H(i) = H(i) - J*M(i, j)*(M(i, j-1) + M(i, j+2));
          elseif j == 3
              H(i) = H(i) - J*M(i, j)*(M(i, j-2) + M(i, j+1));
          else
              H(i) = H(i) - J*M(i, j)*(M(i, j-1) + M(i, j-2));
          end
       else
           if j == 1
              H(i) = H(i) - J*M(i, j)*(M(i, j+1) + M(i, n+1) + M(i, n) + M(i, n^2-n+1));
           elseif j == n
               H(i) = H(i) - J*M(i, j)*(M(i, j-1) + M(i, 1) + M(i, 2*n) + M(i, n^2));
           elseif j == n^2
               H(i) = H(i) - J*M(i, j)*(M(i, j-1) + M(i, n^2-n+1) + M(i, n) + M(i, n^2 - n));
           elseif j == n^2 - n + 1
               H = H - J*M(i, j)*(M(i, j+1) + M(i, n^2) + M(i, n^2-2*n+1) + M(i, 1));
           elseif j < n && j > 1
               H(i) = H(i) - J*M(i, j)*(M(i, j-1) + M(i, j+1) + M(i, n+j) + M(i, n^2-n+j));
           elseif j < n^2 && j > n^2-n+1
               H(i) = H(i) - J*M(i, j)*(M(i, j-1) + M(i, j+1) + M(i, j-n) + M(i, n - (n^2 - j)));
           elseif j ~= 1 && j ~= n^2 - n + 1 && mod(j-1, n) == 0
               H(i) = H(i) - J*M(i, j)*(M(i, j+1) + M(i, j-n) + M(i, j+n) + M(i, j+n-1));
           elseif j ~= n && j ~= n^2 && mod(j, n) == 0
               H(i) = H(i) - J*M(i, j)*(M(i, j-1) + M(i, j-n) + M(i, j-n+1) + M(i, j+n));
           else
               H(i) = H(i) - J*M(i, j)*(M(i, j-1) + M(i, j+1) + M(i, j+n) + M(i, j-n));
           end
       end
    end
    
    H(i) = H(i) - h*sum(M(i, :));
    Z = Z + exp(-H(i)./T);
    He = He + H(i).*exp(-H(i)./T);
    C = C + (H(i)^2.*exp(-H(i)./T));
    for k = 1:n^2
       m = m + M(i, k)*exp(-H(i)./T);
    end
end

C = (C./Z - (He./Z).^2)./(T.^2);
He = He./Z;
m = m./(n^2.*Z);

figure(1)
plot(T, Z)
xlabel('T');
ylabel('Z');
title('Partitioning function');

figure(2)
plot(T, He)
xlabel('T');
ylabel('U');
title('Inner energy');

figure(3)
plot(T, C)
title('Heat capacity');
xlabel('T');
ylabel('C');

figure(4)
plot(T, m)
title('Magnetization');
xlabel('T');
ylabel('m');
hold on