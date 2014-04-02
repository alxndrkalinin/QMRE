% ѕостроение вектора откликов

%     ¬ходные параметры:
% alpha - вектор оцениваемых параметров
% x - регрессоры
% N - количество экспериментов
% J - количество градаций отклика
% z - случайный вектор откликов

%     ¬ыходные параметры:
% z - сгенерированный вектор откликов

% -------------------------------------------------------------------------

function out = CZ(alpha, x, N, J, z)

P = ones(J, 1);
out = zeros(N, 1);
flag = 0;

for i = 1:N
    % ¬ычисление веро€тностей попадани€ в градацию
%     for j = 1:J
%         P(j) = 1;
%         subsum = 1;
%         for k = 1:J-1
%             subsum = subsum + exp(x(i, :) * alpha(k, :));
%         end
%         if j == J
%             pr = 1;
%         else
%             pr =  exp(x(i, :) * alpha(j, :));
%         end
%         P(j) = P(j) * pr / subsum;
%     end

P(1) = exp(x(i, :) * alpha(1, :)) / (1 + exp(x(i, :) * alpha(1, :)) + exp(x(i, :) * alpha(2, :)) + exp(x(i, :) * alpha(3, :)));
P(2) = exp(x(i, :) * alpha(2, :)) / (1 + exp(x(i, :) * alpha(1, :)) + exp(x(i, :) * alpha(2, :)) + exp(x(i, :) * alpha(3, :)));
P(3) = exp(x(i, :) * alpha(3, :)) / (1 + exp(x(i, :) * alpha(1, :)) + exp(x(i, :) * alpha(2, :)) + exp(x(i, :) * alpha(3, :)));
P(4) = 1 / (1 + exp(x(i, :) * alpha(1, :)) + exp(x(i, :) * alpha(2, :)) + exp(x(i, :) * alpha(3, :)));

    % ¬ычисление кумул€тивных веро€тностей
    for j = 1:J-1
        P(j + 1) = P(j) + P(j + 1);
    end

    % ¬ычисление значени€ отклика в зависимости от попадани€ в веро€тности
    for j = 1:J-1
        if z(i) > P(j)
            out(i) = j + 1;
        else
            if flag == 0
                out(i) = j;
                flag = 1;
            end
        end
    end
    flag = 0;
end

end