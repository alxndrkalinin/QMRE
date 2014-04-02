% Вычисление ЛФП

%   Входные параметры:
% alpha - вектор оцениваемых параметров
% x - регрессоры
% z - вектор отклика
% J - размерность вектора отклика
% N - количество экспериментов

%   Выходные параметры:
% out - отрицательное значение функции правдоподобия

% -------------------------------------------------------------------------

function [out, gradL] = CORPSI(alphaj, alpha0, x, z, J, N, lambda, j)

grad = 0;
psisum = 0;

alpha = [alpha0; 0];
alpha(j) = alphaj;
% alpha = alpha0;

% Вычисление Пси
for i = 1:N
    
    psi = 0;
    psi0 = krondelta(j, z(i)) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) / (sumP(x(i, :), alpha, J) ^ ((1 + lambda)));
    sum = 0;
    for m = 1:J
        sum = sum + (P(x(i, :), alpha(m), alpha, J) ^ (1 - lambda));
    end
    psi = psi0 * (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda) * sum * x(i, :);
    psisum = psisum + abs(psi)^2;

end

out = psisum;

% out = -out2;

% % Вычисление градиента ЛПФ
% 
% for s = 1:J-1
%     sum = 0;
%     for i = 1:N
%         ksiz = 0;
%         numsum = 0;
%         denomsum = 0;
%         for j = 1:J
%             ksij = 0;
%             if s == j
%                 ksij = x(i);
%             end
%             numsum = numsum + ksij * exp(x(i,:) * alpha(j));
%         end
%         for j = 1:J
%             denomsum = denomsum + exp(x(i,:) * alpha(j));
%         end
%         if s == z(i)
%             ksiz = 1;
%         end
%         sum = sum + (ksiz - numsum / denomsum);
%     end
%     grad(s) = sum;
% end
gradL = -grad;

end

function [prob] = P(x, alpha_s, alpha, J)
    
%    prob = 1.1 ^ (x * alpha_s) / sumP(x, alpha, J);
    prob = exp(x * alpha_s) / sumP(x, alpha, J);
    
end

function [sumprob] = sumP(x, alpha, J)

    sumprob = 0;
    for m = 1:J
        sumprob = sumprob + exp(x * alpha(m));
%        sumprob = sumprob + 1.1 ^(x * alpha(m));
    end

end

function [out] = krondelta(a, b)

    out = 0;
    if a == b
        out = 1;
    end

end