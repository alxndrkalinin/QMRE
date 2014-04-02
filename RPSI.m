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

function [out, gradL] = RPSI(alphaj, alpha0, x, z, J, N, lambda, j)

grad = zeros(1, J);
psisum = 0;
gradL = 0;

% alpha = [alpha0; 0];
alpha = zeros()
alpha(j) = alphaj;
test = J;
% alpha = alpha0;

% Вычисление Пси
for i = 1:N
    
    psi0 = krondelta(j, z(i)) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) / (expsum(x(i, :), alpha, J) ^ (1 + lambda));
    delta = deltaP(x(i, :), alpha, J, lambda);
    psi = psi0 * (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda) * delta * x(i, :);
    psisum = psisum + abs(psi)^2.2;
    
%     for l = 1:J-1
%         grad(l) = grad(l) + ((- ((1 + lambda) * (derP(x(i, :), alpha, J, j, l) ^ lambda) * (sumP(x(i, :), i, z, alpha, J, (1 + lambda))) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) * ((1 + lambda) * dersumP(x(i, :), alpha, J, j, l, (lambda)))) / (sumP(x(i, :), i, z, alpha, J, (1 + lambda)) ^ 2)) * (P(x(i, :), alpha(j), alpha, J) ^ lambda) + (krondelta(j, z(i)) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) / (sumP(x(i, :), i, z, alpha, J, (1 + lambda)))) * (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda)) * (deltaP(x(i, :), alpha, J, lambda)) + ((krondelta(j, z(i)) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) / (sumP(x(i, :), i, z, alpha, J, (1 + lambda)))) * (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda)) * (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda) * (sumP(x(i, :), i, z, alpha, J, (1 - lambda)));
%     end
    l = j;
    gradL = gradL + ((- ((1 + lambda) * (derP(x(i, :), alpha, J, j, l) ^ lambda) * (sumP(x(i, :), i, z, alpha, J, (1 + lambda))) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) * ((1 + lambda) * dersumP(x(i, :), alpha, J, j, l, (lambda)))) / (sumP(x(i, :), i, z, alpha, J, (1 + lambda)) ^ 2)) * (P(x(i, :), alpha(j), alpha, J) ^ lambda) + (krondelta(j, z(i)) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) / (sumP(x(i, :), i, z, alpha, J, (1 + lambda)))) * (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda)) * (deltaP(x(i, :), alpha, J, lambda)) + ((krondelta(j, z(i)) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) / (sumP(x(i, :), i, z, alpha, J, (1 + lambda)))) * (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda)) * (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda) * (sumP(x(i, :), i, z, alpha, J, (1 - lambda)));
end

out = psisum;
% gradL = grad(j);

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
% gradL = -grad;

end

function [prob] = P(x, alpha_s, alpha, J)
    
   prob = (1.1 ^ (x * alpha_s)) / expsum(x, alpha, J);
%     prob = exp(x * alpha_s) / expsum(x, alpha, J);
    
end

function [sumprob] = expsum(x, alpha, J)

    sumprob = 0;
    for m = 1:J-1
%        sumprob = sumprob + exp(x * alpha(m));
        sumprob = sumprob + 1.1 ^(x * alpha(m));
    end

end

function [sumprob] = sumP(x, i, z, alpha, J, power)

    sumprob = 0;
    for m = 1:J-1
        sumprob = sumprob + P(x, alpha(z(i)), alpha, J) ^ power;
    end

end

function [out] = krondelta(a, b)

    out = 0;
    if a == b
        out = 1;
    end

end

function [out] = deltaP(x, alpha, J, lambda)
    out = 0;
    for m = 1:J-1
        out = out + (P(x, alpha(m), alpha, J) ^ (1 - lambda));
    end
%     out = 1;
end

% s - индекс оцениваемой компоненты альфа, k - по которой берем производную
function [derprob] = derP(x, alpha, J, s, k)
    
    % сумма, под которой крондельта
    deltedsum = 0;
    for m = 1:J-1
%         deltedsum = deltedsum + krondelta(alpha(m), alpha(k)) * exp(x * alpha(m));
        deltedsum = deltedsum + krondelta(alpha(m), alpha(k)) * (1.1 ^ (x * alpha(m)));
    end
%     derprob = x * exp(x * alpha(s)) * ((krondelta(alpha(s), alpha(k)) * (1 + expsum(x, alpha, J)) - deltedsum)) / ((1 + expsum(x, alpha, J)) ^ 2);
    derprob = x * (1.1 ^ (x * alpha(s))) * ((krondelta(alpha(s), alpha(k)) * (1 + expsum(x, alpha, J)) - deltedsum)) / ((1 + expsum(x, alpha, J)) ^ 2);
    
end

function [deltedsum] = dersumP(x, alpha, J, s, k, power)
    
    % сумма, под которой крондельта
    deltedsum = 0;
    for m = 1:J-1
        deltedsum = deltedsum + derP(x, alpha, J, s, k) ^ power;
    end    
end