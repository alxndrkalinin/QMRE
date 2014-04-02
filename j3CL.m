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

function [out, gradL] = j3CL(alpha0, x, z, J, N)

out2 = 0;
grad = 0;
a=[0;0;0];

alpha = [alpha0, a];
% alpha = alpha0;

% Вычисление ЛПФ
for i = 1:N
    out1 = 0;
    for j = 1:J
        xx = [1, x(i,:), x(i,:)^2];
        out1 = out1 + exp(sum(xx * alpha(:,j)));
    end
    out2 = out2 + sum(xx * alpha(:,z(i))) - log(out1);
end

out = -out2;

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