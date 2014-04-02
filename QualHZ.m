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

function QQ = QualHZ(alpha1, alpha2, x, N, hJ, z, hJprod, hJshift, hnz)

P1 = zeros(hnz, hJprod);
P2 = zeros(hnz, hJprod);
out1 = zeros(N, hnz);
out2 = zeros(N, hnz);
QQ = 0;
flag1 = 0;
flag2 = 0;

for i = 1:N
%     Q = 0;
%     % Вычисление вероятностей для первого отклика
%     for j = 1:hJ(1)
%             P1(1, j) = 1;
%             P2(1, j) = 1;
%             subsum1 = 1;
%             subsum2 = 1;
%             for k = 1:hJ(1)-1
%                 subsum1 = subsum1 + exp(x(i, :) * alpha1(k, :));
%                 subsum2 = subsum2 + exp(x(i, :) * alpha2(k, :));
%             end
%             if j == hJ(1)
%                 pr1 = 1;
%                 pr2 = 1;
%             else
%                 pr1 =  exp(x(i, :) * alpha1(j, :));
%                 pr2 =  exp(x(i, :) * alpha2(j, :));
%             end
%             P1(1, j) = P1(1, j) * pr1 / subsum1;
%             P2(1, j) = P2(1, j) * pr2 / subsum2;
% %             Q = Q + (P1(1, j) - P2(1, j))^2;
%     end
%          % Вычисление кумулятивных вероятностей
%          Pcum1 = P1;
%          Pcum2 = P2;
%     for j = 1:hJ(1)-1
%         Pcum1(1, j + 1) = Pcum1(1, j) + Pcum1(1, j + 1);
%         Pcum2(1, j + 1) = Pcum2(1, j) + Pcum2(1, j + 1);
%     end
% 
%     % Вычисление значения отклика в зависимости от попадания в вероятности
%     for j = 1:hJ(1)-1
%         if z(i) > Pcum1(1, j)
%             out1(i, 1) = j + 1;
%         else
%             if flag1 == 0
%                 out1(i, 1) = j;
%                 flag1 = 1;
%             end
%         end
%         if z(i) > Pcum2(1, j)
%             out2(i, 1) = j + 1;
%         else
%             if flag2 == 0
%                 out2(i, 1) = j;
%                 flag2 = 1;
%             end
%         end
%     end
%     flag1 = 0;
%     flag2 = 0;
%     QQ = QQ + Q;
% % -------------------------------------------------------------------------
%     Q = 0;
%     for k = 2:hnz                               % цикл по откликам
%         Ptmp1 = ones(1, hJ(k));                  % вероятности отклика
%         Ptmp2 = ones(1, hJ(k));                  % вероятности отклика
%         outshift1 = out1(i, k - 1);               % сдвиг по парамерам
%         outshift2 = out2(i, k - 1);
%         for kk = 3:k
%             outshift1 = outshift1 + (out1(i, kk - 1) - 1) * hJprod(kk - 2);
%             outshift2 = outshift2 + (out2(i, kk - 1) - 1) * hJprod(kk - 2);
%         end
%         for j = 1:hJ(k)                       % цикл по градациям
%             % Вычисление своей вероятности
%             subsum1 = 1;
%             subsum2 = 1;
%             for s = 1:hJ(k)-1
%                 subsum1 = subsum1 + exp(x(i, :) * alpha1(hJshift(k) + outshift1 + hJprod(k - 1) * (s - 1), :));
%                 subsum2 = subsum2 + exp(x(i, :) * alpha2(hJshift(k) + outshift2 + hJprod(k - 1) * (s - 1), :));
%             end
%             if j == hJ(k)
%                 pr1 = 1;
%                 pr2 = 1;
%             else
%                 pr1 =  exp(x(i, :) * alpha1(hJshift(k) + outshift1 + hJprod(k - 1) * (j - 1), :));
%                 pr2 =  exp(x(i, :) * alpha2(hJshift(k) + outshift2 + hJprod(k - 1) * (j - 1), :));
%             end
%             Ptmp1(j) = Ptmp1(j) * pr1 / subsum1;
%             Ptmp2(j) = Ptmp2(j) * pr2 / subsum2;  
%         end
%         P1(k, :) = Ptmp1;
%         P2(k, :) = Ptmp2;
%         
%         for jj = 1:hJ(k)
%             Pprod1(k, (jj-1)*hJ(k-1)+1:jj*hJ(k-1)) = P1(k, jj) * P1(k-1, :);
%             Pprod2(k, (jj-1)*hJ(k-1)+1:jj*hJ(k-1)) = P2(k, jj) * P2(k-1, :);
%         end
%         
%          % Вычисление кумулятивных вероятностей
%          Pcum1(k, :) = P1(k, :);
%          Pcum2(k, :) = P2(k, :);
%         for j = 1:hJ(k)-1
%             Pcum1(k, j + 1) = Pcum1(k, j) + Pcum1(k, j + 1);
%             Pcum2(k, j + 1) = Pcum2(k, j) + Pcum2(k, j + 1);
%         end
% 
%         % Вычисление значения отклика в зависимости от попадания в вероятности
%         for j = 1:hJ(k)-1
%             if z(i) > Pcum1(k, j)
%                 out1(i, k) = j + 1;
%             else
%                 if flag1 == 0
%                     out1(i, k) = j;
%                     flag1 = 1;
%                 end
%             end
%             if z(i) > Pcum2(k, j)
%                 out2(i, k) = j + 1;
%             else
%                 if flag2 == 0
%                     out2(i, k) = j;
%                     flag2 = 1;
%                 end
%             end
%         end
%     end
    
    Q=0;
    
A1 = exp(x(i, :) * alpha1(1)) / (1 + exp(x(i, :) * alpha1(1))) * exp(x(i, :) * alpha1(2)) / (1 + exp(x(i, :) * alpha1(2)));
A2 = exp(x(i, :) * alpha1(1)) / (1 + exp(x(i, :) * alpha1(1))) * 1 / (1 + exp(x(i, :) * alpha1(2)));
A3 = 1 / (1 + exp(x(i, :) * alpha1(1))) * exp(x(i, :) * alpha1(3)) / (1 + exp(x(i, :) * alpha1(3)));
A4 = 1 / (1 + exp(x(i, :) * alpha1(1))) * 1 / (1 + exp(x(i, :) * alpha1(3)));

Pprod1(hnz,1)=A1;
Pprod1(hnz,2)=A2;
Pprod1(hnz,3)=A3;
Pprod1(hnz,4)=A4;

B1 = exp(x(i, :) * alpha2(1)) / (1 + exp(x(i, :) * alpha2(1))) * exp(x(i, :) * alpha2(2)) / (1 + exp(x(i, :) * alpha2(2)));
B2 = exp(x(i, :) * alpha2(1)) / (1 + exp(x(i, :) * alpha2(1))) * 1 / (1 + exp(x(i, :) * alpha2(2)));
B3 = 1 / (1 + exp(x(i, :) * alpha2(1))) * exp(x(i, :) * alpha2(3)) / (1 + exp(x(i, :) * alpha2(3)));
B4 = 1 / (1 + exp(x(i, :) * alpha2(1))) * 1 / (1 + exp(x(i, :) * alpha2(3)));

Pprod2(hnz,1)=B1;
Pprod2(hnz,2)=B2;
Pprod2(hnz,3)=B3;
Pprod2(hnz,4)=B4;
    
    for jj = 1:hJprod(hnz)
            Q = Q + (Pprod1(hnz,jj) - Pprod2(hnz,jj))^2;
    end
    
    flag1 = 0;
    flag2 = 0;
    QQ = QQ + Q;
end
% QQ = QQ / (N * hJprod(hnz));
QQ = QQ / N;
end