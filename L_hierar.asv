% Вычисление ЛФП

%   Входные параметры:
% alpha - вектор оцениваемых параметров
% x - регрессоры
% z - вектор отклика
% J - размерность вектора отклика
% N - количество экспериментов

%   Выходные параметры:
% out - отрицательное значение функции правдоподобия

function out = L_hierar(alpha0, x, z, J, N)

out2 = 0;

%alpha = [alpha0(end-1); 0];
alpha = alpha0;

% Вычисление ЛПФ
for i = 1:N
    out1 = 0;
    for j = 1:J
        out1 = out1 + exp(x(i,:) * alpha(j));
    end
    out2 = out2 + x(i,:) * alpha(z(i)) - log(out1);
end

out = -out2;