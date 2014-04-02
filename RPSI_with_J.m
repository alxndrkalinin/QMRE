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

function [out, globalgrad] = RPSI_with_J(alpha0, x, z, J, N, lambda, j)

globalgrad = zeros(1, J-1);
out = 0;
alpha = [alpha0; 0];

% alpha = [alpha0; 0];
% alpha = zeros()
% alpha(j) = alphaj;
% alpha = alpha0;

% Вычисление Пси
    for j = 1:J-1       %  по оцениваемой компоненте альфа
        psisum = 0;
        localgrad = zeros(1, J-1);
        
        for l = 1:J-1       %  по какой дифференцируем

            gradPsi = 0;
            for i = 1:N

                if l == 1
                    psi0 = krondelta(j, z(i)) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) / (expsum(x(i, :), alpha, J) ^ (1 + lambda));
                    delta = deltaP(x(i, :), alpha, J, lambda);
                    psi = psi0 * (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda) * delta * x(i, :);
                    psisum = psisum + psi;
                end

        %         for l = 1:J-1

        %             grad(l) = grad(l) + ((- ((1 + lambda) * (derP(x(i, :), alpha, J, j, l) ^ lambda) * (sumP(x(i, :), i, z, alpha, J, (1 + lambda))) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) * ((1 + lambda) * dersumP(x(i, :), alpha, J, j, l, (lambda)))) / (sumP(x(i, :), i, z, alpha, J, (1 + lambda)) ^ 2)) * (P(x(i, :), alpha(j), alpha, J) ^ lambda) + (krondelta(j, z(i)) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) / (sumP(x(i, :), i, z, alpha, J, (1 + lambda)))) * (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda)) * (deltaP(x(i, :), alpha, J, lambda)) + ((krondelta(j, z(i)) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) / (sumP(x(i, :), i, z, alpha, J, (1 + lambda)))) * (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda)) * (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda) * (sumP(x(i, :), i, z, alpha, J, (1 - lambda)));
        %         end
        %         l = j;
        %         gradL = gradL + ((- ((1 + lambda) * (derP(x(i, :), alpha, J, j, l) ^ lambda) * (sumP(x(i, :), i, z, alpha, J, (1 + lambda))) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) * ((1 + lambda) * dersumP(x(i, :), alpha, J, j, l, (lambda)))) / (sumP(x(i, :), i, z, alpha, J, (1 + lambda)) ^ 2)) * (P(x(i, :), alpha(j), alpha, J) ^ lambda) + (krondelta(j, z(i)) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) / (sumP(x(i, :), i, z, alpha, J, (1 + lambda)))) * (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda)) * (deltaP(x(i, :), alpha, J, lambda)) + ((krondelta(j, z(i)) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) / (sumP(x(i, :), i, z, alpha, J, (1 + lambda)))) * (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda)) * (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda) * (sumP(x(i, :), i, z, alpha, J, (1 - lambda)));
    %             gradL = gradL + ((1 + lambda) * ((P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) * (dermixsumP(x, alpha, J, j, l, 1, lambda)) - (derP(x, alpha, J, l, k)) * () * ()) / ()) * () * () + () * ();
                gradPsi = gradPsi + (((1 + lambda) * ((P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) * (dermixsumP(x(i, :), alpha, J, l, 1, lambda)) - ...
                    (derP(x(i, :), alpha, J, j, l)) * (P(x(i, :), alpha(j), alpha, J) ^ lambda) * (sumP(x(i, :), i, z, alpha, J, (1 + lambda)))) / ...
                    ((sumP(x(i, :), i, z, alpha, J, (1 + lambda))) ^ 2)) * (P(x(i, :), alpha(z(i)), alpha, J)) * (sumP(x(i, :), i, z, alpha, J, (1 - lambda))) + ...
                    (krondelta(j, z(i)) - (P(x(i, :), alpha(j), alpha, J) ^ (1 + lambda)) / (sumP(x(i, :), i, z, alpha, J, (1 + lambda)))) * (lambda * derP(x(i, :), alpha, J, z(i), l) * ...
                    (P(x(i, :), alpha(z(i)), alpha, J) ^ (lambda - 1)) * (sumP(x(i, :), i, z, alpha, J, (1 - lambda))) + (1 - lambda) * ...
                    (P(x(i, :), alpha(z(i)), alpha, J) ^ lambda) * (dermixsumP(x(i, :), alpha, J, l, 1, -lambda))));
            end

            localgrad(l) = localgrad(l) + psisum * gradPsi;
        end
        
        out = out + abs(psisum) ^ 2;
        globalgrad = globalgrad + localgrad;
    end

end

function [prob] = P(x, alpha_s, alpha, J)
    
%    prob = (1.05 ^ (x * alpha_s)) / expsum(x, alpha, J);
    prob = exp(x * alpha_s) / expsum(x, alpha, J);
    
end

function [sumprob] = expsum(x, alpha, J)

    sumprob = 0;
    for m = 1:J-1
       sumprob = sumprob + exp(x * alpha(m));
%         sumprob = sumprob + 1.05 ^(x * alpha(m));
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
%     out = 0;
%     for m = 1:J-1
%         out = out + (P(x, alpha(m), alpha, J) ^ (1 - lambda));
%     end
    out = 1;
end

% s - индекс оцениваемой компоненты альфа, k - по которой берем производную
function [derprob] = derP(x, alpha, J, s, k)
    
    % сумма, под которой крондельта
%     deltedsum = 0;
%     for m = 1:J-1
%         deltedsum = deltedsum + krondelta(alpha(m), alpha(k)) * exp(x * alpha(m));
%         deltedsum = deltedsum + krondelta(alpha(m), alpha(k)) * (1.1 ^ (x * alpha(m)));
%     end
%     derprob = x * exp(x * alpha(s)) * ((krondelta(alpha(s), alpha(k)) * (1 + expsum(x, alpha, J)) - deltedsum)) / ((1 + expsum(x, alpha, J)) ^ 2);
%     derprob = x * (1.1 ^ (x * alpha(s))) * ((krondelta(alpha(s), alpha(k)) * (1 + expsum(x, alpha, J)) - deltedsum)) / ((1 + expsum(x, alpha, J)) ^ 2);
    derprob = P(x, alpha(s), alpha, J) * (krondelta(alpha(s), alpha(k)) - P(x, alpha(k), alpha, J)) * x;
end

% индекс k - по которой берем производную
function [deltedsum] = dermixsumP(x, alpha, J, k, derpower, power)
    
    % сумма, под которой крондельта
    deltedsum = 0;
    for m = 1:J-1
        deltedsum = deltedsum + (derP(x, alpha, J, m, k) ^ derpower) *  (P(x, alpha(m), alpha, J) ^ power);
    end    
end