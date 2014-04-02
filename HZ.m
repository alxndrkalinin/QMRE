% Построение вектора откликов

%     Входные параметры:
% alpha - вектор оцениваемых параметров
% x - регрессоры
% N - количество экспериментов
% J - количество градаций отклика
% z - случайный вектор откликов
% P - массив вероятностей размером hnz * hJprod

%     Выходные параметры:
% z - сгенерированный вектор откликов

% -------------------------------------------------------------------------

function out = HZ(alpha, x, N, hJ, z, hJprod, hJshift, hnz)

P = zeros(hnz, hJprod);
out = zeros(N, hnz);
flag = 0;

for i = 1:N

%     % Вычисление вероятностей для первого отклика
%     for j = 1:hJ(1)
%             P(1, j) = 1;
%             subsum = 1;
%             for k = 1:hJ(1)-1
%                 subsum = subsum + exp(x(i, :) * alpha(k, :));
%             end
%             if j == hJ(1)
%                 pr = 1;
%             else
%                 pr =  exp(x(i, :) * alpha(j, :));
%             end
%             P(1, j) = P(1, j) * pr / subsum;
%     end
%     % Вычисление кумулятивных вероятностей
%     Pcum = P;
%     for j = 1:hJ(1)-1
%         Pcum(1, j + 1) = Pcum(1, j) + Pcum(1, j + 1);
%     end
% 
%     % Вычисление значения отклика в зависимости от попадания в вероятности
%     for j = 1:hJ(1)-1
%         if z(i) > Pcum(1, j)
%             out(i, 1) = j + 1;
%         else
%             if flag == 0
%                 out(i, 1) = j;
%                 flag = 1;
%             end
%         end
%     end
%     flag = 0;
% % -------------------------------------------------------------------------
%     for k = 2:hnz                               % цикл по откликам
%         Ptmp = ones(1, hJ(k));                  % вероятности отклика
%         outshift = out(i, k - 1);               % сдвиг по парамерам
%         for kk = 3:k
%             outshift = outshift + (out(i, kk - 1) - 1) * hJprod(kk - 2);
%         end
%         for j = 1:hJ(k)                       % цикл по градациям
%             % Вычисление своей вероятности
%             subsum = 1;
%             for s = 1:hJ(k)-1
%                 subsum = subsum + exp(x(i, :) * alpha(hJshift(k) + outshift + hJprod(k - 1) * (s - 1), :));
%             end
%             if j == hJ(k)
%                 pr = 1;
%             else
%                 pr =  exp(x(i, :) * alpha(hJshift(k) + outshift + hJprod(k - 1) * (j - 1), :));
%             end
%             Ptmp(j) = Ptmp(j) * pr / subsum;
%         end
%         
%         P(k, :) = Ptmp;
        
A1 = exp(x(i, :) * alpha(1)) / (1 + exp(x(i, :) * alpha(1))) * exp(x(i, :) * alpha(2)) / (1 + exp(x(i, :) * alpha(2)));
A2 = exp(x(i, :) * alpha(1)) / (1 + exp(x(i, :) * alpha(1))) * 1 / (1 + exp(x(i, :) * alpha(2)));
A3 = 1 / (1 + exp(x(i, :) * alpha(1))) * exp(x(i, :) * alpha(3)) / (1 + exp(x(i, :) * alpha(3)));
A4 = 1 / (1 + exp(x(i, :) * alpha(1))) * 1 / (1 + exp(x(i, :) * alpha(3)));

Pprod1(hnz,1)=A1;
Pprod1(hnz,2)=A2;
Pprod1(hnz,3)=A3;
Pprod1(hnz,4)=A4;

    % Вычисление кумулятивных вероятностей
    for j = 1:hJprod(hnz)-1
        Pprod1(hnz, j + 1) = Pprod1(hnz, j) + Pprod1(hnz, j + 1);
    end

    % Вычисление значения отклика в зависимости от попадания в вероятности
    for j = 1:hJprod(hnz)-1
        if z(i) > Pprod1(hnz, j)
            out(i, 1) = j + 1;
        else
            if flag == 0
                out(i, 1) = j;
                flag = 1;
            end
        end
    end
     flag = 0;
end

        
%         % Вычисление кумулятивных вероятностей
%         Pcum(k, :) = P(k, :);
%         for j = 1:hJ(k)-1
%             Pcum(k, j + 1) = Pcum(k, j) + Pcum(k, j + 1);
%         end
% 
%         % Вычисление значения отклика в зависимости от попадания в вероятности
%         for j = 1:hJ(k)-1
%             if z(i) > Pcum(k, j)
%                 out(i, k) = j + 1;
%             else
%                 if flag == 0
%                     out(i, k) = j;
%                     flag = 1;
%                 end
%             end
%         end
%     end
    
   
    out = out(:,1);
end



% end