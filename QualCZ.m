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

function out = QualCZ(alpha1, alpha2, x, N, J, z)

P1 = ones(J, 1);
P2 = ones(J, 1);
out = 0;
Q = 0;

for i = 1:N
    % ¬ычисление веро€тностей попадани€ в градацию
    Q = 0;
%     for j = 1:J
%         P1(j) = 1;
%         P2(j) = 1;
%         subsum1 = 1;
%         subsum2 = 1;
%         for k = 1:J-1
%             subsum1 = subsum1 + exp(x(i, :) * alpha1(k, :));
%             subsum2 = subsum2 + exp(x(i, :) * alpha2(k, :));
%         end
%         if j == J
%             pr1 = 1;
%             pr2 = 1;
%         else
%             pr1 =  exp(x(i, :) * alpha1(j, :));
%             pr2 =  exp(x(i, :) * alpha2(j, :));
%         end
%         P1(j) = P1(j) * pr1 / subsum1;
%         P2(j) = P2(j) * pr2 / subsum2;
%         Q = Q + (P1(j) - P2(j))^2;
%     end

    P1(1) = exp(x(i, :) * alpha1(1, :)) / (1 + exp(x(i, :) * alpha1(1, :)) + exp(x(i, :) * alpha1(2, :)) + exp(x(i, :) * alpha1(3, :)));
    P1(2) = exp(x(i, :) * alpha1(2, :)) / (1 + exp(x(i, :) * alpha1(1, :)) + exp(x(i, :) * alpha1(2, :)) + exp(x(i, :) * alpha1(3, :)));
    P1(3) = exp(x(i, :) * alpha1(3, :)) / (1 + exp(x(i, :) * alpha1(1, :)) + exp(x(i, :) * alpha1(2, :)) + exp(x(i, :) * alpha1(3, :)));
    P1(4) = 1 / (1 + exp(x(i, :) * alpha1(1, :)) + exp(x(i, :) * alpha1(2, :)) + exp(x(i, :) * alpha1(3, :)));

    P2(1) = exp(x(i, :) * alpha2(1, :)) / (1 + exp(x(i, :) * alpha2(1, :)) + exp(x(i, :) * alpha2(2, :)) + exp(x(i, :) * alpha2(3, :)));
    P2(2) = exp(x(i, :) * alpha2(2, :)) / (1 + exp(x(i, :) * alpha2(1, :)) + exp(x(i, :) * alpha2(2, :)) + exp(x(i, :) * alpha2(3, :)));
    P2(3) = exp(x(i, :) * alpha2(3, :)) / (1 + exp(x(i, :) * alpha2(1, :)) + exp(x(i, :) * alpha2(2, :)) + exp(x(i, :) * alpha2(3, :)));
    P2(4) = 1 / (1 + exp(x(i, :) * alpha2(1, :)) + exp(x(i, :) * alpha2(2, :)) + exp(x(i, :) * alpha2(3, :)));

    for j =1:J
    Q = Q + (P1(j) - P2(j))^2;
    end
    out = out + Q;
end

% out = out / (N * J);
out = out / N;

end