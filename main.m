% Главный модуль программы

% Оценивание параметров моделей

% -------------------------------------------------------------------------

% Чистка
clear
clc
format short g

M = 1;        % число экспериментов Монте-Карло

% Инициализация данных о моделях
[N, Xdim, XWdim, xmaxpow, cJ, hnz, hJ] = params();
N, Xdim, XWdim, xmaxpow, cJ

x(1,1) = 1;
x(1,2) = -10;
step = 20 / N;
for i = 2:N
    x(i,1) = 1;
    x(i,2) = x(i-1,2) + step;
end
% x(:,3) = x(:,2);

x = x(:,2);

lambda = 1;
eps = 0.1

% x = -ones(N, 1);

alpha1 = zeros(cJ-1, M);
alphar = zeros(cJ-1, M);
alphaco = zeros(cJ-1, M);
% alpha2 = zeros(cJ-1, M);
% alpha3 = zeros(cJ-1, M);
% alpha4 = zeros(cJ-1, M);

Lval1 = zeros(1, M);
% Lval2 = zeros(1, M);
% Lval3 = zeros(1, M);
% Lval4 = zeros(1, M);

% qual1 = zeros(N);
% qual2 = zeros(N);
% qual3 = zeros(N);
% qual4 = zeros(N);

% creal = [-0.5; 0.5; 1.5; 1.6; -0.5; 1.3; -0.7; 1];
% hreal = [-0.5; 0.5; 1.5; 1.6; -0.5; 1.3; -0.7; 1];\

% creal = [-0.5; 0.5; 1];
% rreal = [-0.5; 0.5; 1];

% coreal = [-0.5; 0.5; 1];
% hreal = [-0.5; -1; 1];

creal = [-8, -5; 2, 4; 1, 1]
rreal = creal; coreal = creal;

% creal = [1; 4; 3];
% hreal = [1; -3; 3];

% creal = [2; 0.2; 0.1];
% hreal = [-1; 0.1; 2];

for i = 1:M

% -------------------------------------------------------------------------
% КЛАССИЧЕСКАЯ МОДЕЛЬ

z = rand(N, 1);

% Истинный вектор параметров классической модели
% calpha = [1; 4; 3]
% calpha = [-0.5; 0.5; 1]
calpha = [-8, -5; 2, 4; 1, 1]
% calpha = [0.5; -1; 1; -0.5; 3; -2; 2; -1];
% calpha = [2; 0.2; 0.1]

% Построение ветора откликов классической модели
cz = j3CZ(calpha, x, N, cJ, z(:, 1), lambda, eps);
% Вычисление функции правдоподобия от истинных параметров классической
% модели
% cLreal = - CL(calpha, x, cz, cJ, N)

fid = fopen('z.txt', 'wt');
fprintf(fid, '%f\n', cz);
fclose(fid);

fid = fopen('x.txt', 'wt');
fprintf(fid, '%f\n', x);
fclose(fid);

% Начальное приближение параметров и ЛФП от начального приближения
% классической модели
% calpha0 = [0;0;0];
calpha0 = [0, 0; 0, 0; 0, 0;]
% calpha0 = [0;0;0;0;0;0;0;0];
% cL0 = - CL(calpha0, x, cz, cJ, N)

% Поиск максимума функции правдоподобия классической модели
i
options = optimset('LargeScale', 'off', 'GradObj', 'off', 'Display', 'iter', 'MaxIter', 120, 'MaxFunEvals', 1000);
% [calpha, cFVAL, cEXITFLAG, coutput, cgrad, chessian] = fminunc('CL', calpha0, options, x, cz, cJ, N);
[calpha, cFVAL] = fminunc('j3CL', calpha0, options, x, cz, cJ, N)
% quality_calpha = QualCZ(calpha, creal, x, N, cJ, z(:, 1))
% calpha, cFVAL, cEXITFLAG, coutput, cgrad, chessian

% -------------------------------------------------------------------------
% РОБАСТНАЯ МОДЕЛЬ

z = rand(N, 1);
% Истинный вектор параметров классической модели
% calpha = [1; 4; 3]
% ralpha = [-0.5; 0.5; 1]
% coalpha = ralpha
% calpha = [0.5; -1; 1; -0.5; 3; -2; 2; -1];
% calpha = [2; 0.2; 0.1]

rJ = cJ;
coJ = cJ;

% lb = calpha - 2;
% ub = calpha + 2;
% lb = [-10, -6; 0, 0; 0, 0];
% ub = [0, 0; 5, 6; 3, 3];
lb = [-10, -10; -10, -10; -10, -10];
ub = [10, 10; 10, 10; 10, 10];

% Построение ветора откликов классической модели
% rz = RZ(ralpha, x, N, rJ, z(:, 1));
rz = cz;
coz = cz;
% Вычисление функции правдоподобия от истинных параметров классической
% модели
% cLreal = - CL(calpha, x, cz, cJ, N)

% Начальное приближение параметров и ЛФП от начального приближения
% классической модели
%  ralpha0 = [1.5; 1.5; 1.5];
%  coalpha0 = ralpha0;
%  coalpha0 = [0;0;0];
%  ralpha0 = calpha;

% calpha0 = [0;0;0;0;0;0;0;0];
% cL0 = - CL(calpha0, x, cz, cJ, N)

% Поиск радикальных оценок
i
% options = optimset('LargeScale', 'off', 'GradObj','off', 'Display', 'iter', 'MaxIter', 10);
% [calpha, cFVAL, cEXITFLAG, coutput, cgrad, chessian] = fminunc('CL', calpha0, options, x, cz, cJ, N);

% for j = 1:rJ-1
%     options = optimset('LargeScale', 'off', 'GradObj','off', 'Display', 'iter', 'MaxIter', 10);
% %    [ralpha_j, cFVAL] = fminunc('revRPSI', ralpha0(j), options, ralpha0, x, cz, cJ, N, lambda, j)
%     [ralpha_j] = fmincon('revRPSI', ralpha0(j), [], [], [], [], lb, ub, [], options, ralpha0, x, cz, cJ, N, lambda, j);
%     ralpha0(j) = ralpha_j;
%     options = optimset('LargeScale', 'on', 'GradObj','on', 'Display', 'iter', 'MaxIter', 10);
% %    [ralpha_j, cFVAL] = fminunc('RPSI', ralpha0(j), options, ralpha0, x, cz, cJ, N, lambda, j)
%     [ralpha_j, cFVAL] = fmincon('RPSI', ralpha0(j), [], [], [], [], lb, ub, [], options, ralpha0, x, cz, cJ, N, lambda, j)
%     ralpha(j) = ralpha_j;
% end
% 
% quality_ralpha = QualCZ(ralpha, rreal, x, N, cJ, z(:, 1))
% calpha, cFVAL, cEXITFLAG, coutput, cgrad, chessian

% Поиск условно-оптимальных оценок
i
options = optimset('LargeScale', 'on', 'GradObj','on', 'Display', 'iter', 'MaxIter', 100);
 coalpha0 = calpha;
% [calpha, cFVAL, cEXITFLAG, coutput, cgrad, chessian] = fminunc('CL', calpha0, options, x, cz, cJ, N);
% for j = 1:rJ-1
%    [ralpha_j, cFVAL] = fminunc('revRPSI', ralpha0(j), options, ralpha0, x, cz, cJ, N, lambda, j)
%    [ralpha_j] = fmincon('revRPSI', ralpha0(j), [], [], [], [], lb, ub, [], options, ralpha0, x, cz, cJ, N, lambda, j);
%    ralpha0(j) = ralpha_j;
%     [coalpha, cFVAL] = fminunc('RPSI_with_J', coalpha0, options, x, coz,
%     coJ, N, lambda, j)
    
%     for j = 1:cJ-1
%         calpha_lowb = calpha0(:,j);
%         calpha_est = coalpha0(:,j);
%         [int_calpha, cFVAL] = fmincon('Int_with_J', calpha_est, [], [], [], [], lb, ub, [], options, calpha_lowb, x, coz, coJ, N, lambda, j)
%         coalpha0(:, j) = int_calpha;        
%     end

   [coalpha, cFVAL] = fmincon('j3RPSI_with_J', coalpha0, [], [], [], [], lb, ub, [], options, x, coz, coJ, N, lambda, j)
   
%    options = optimset('LargeScale', 'off', 'GradObj','off', 'Display', 'iter', 'MaxIter', 100);
%    coalpha0 = coalpha;
%    [coalpha, cFVAL] = fmincon('j3RPSI_with_J', coalpha0, [], [], [], [], lb, ub, [], options, x, coz, coJ, N, lambda, j)
%     coalpha(j) = coalpha_j;
% end

% quality_coalpha = QualCZ(coalpha, coreal, x, N, cJ, z(:, 1))

% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% % ИЕРАРХИЧЕСКАЯ МОДЕЛЬ
% 
% z = rand(N, hnz);
% 
% % Истинный вектор параметров иерахической модели
% hJprod = hJ;                    % произведение количества всех градаций
% halpha_dim = 0;                 % размерность вектора параметров
% hJshift = hJ;        % сдвиг для вычисления позиции параметра
% for k = 1:hnz
%     if k == 1
%         halpha_dim = hJ(1) - 1;
%         hJshift(k) = 0;
%         hJshift(k + 1) = hJ(1) - 1;
%     else
%         halpha_dim = halpha_dim + (hJ(k) - 1) * hJprod(k - 1);
%         hJprod(k) = hJprod(k) * hJprod(k - 1);
%         hJshift(k + 1) = hJshift(k) + hJprod(k - 1) * (hJ(k) - 1);
%     end
% end
% z = rand(N, hnz);
% 
% % halpha = ones(halpha_dim, Xdim);
% 
% % halpha = [1; -3; 3]
% halpha = [-0.5; -1; 1];
% % halpha = [0.5; -1; 1; -0.5; 3; -2; 2; -1]
% % halpha = [-1; 0.1; 2]
% 
% hJprod
% 
% % Построение ветора откликов иерахической модели
% cz_pr = HZ(halpha, x, N, hJ, z, hJprod, hJshift, hnz);
% 
% % -------------------------------------------------------------------------
% % Перестроение векторов откликов для иерархич модели
% hz = cz_pr;
% for j = 1:N
%     shift = 0;
%     for zz = 1:hnz - 1
%        hz(j, zz) = ceil((cz_pr(j) - shift) / hJprod(hnz - zz));
%        shift = shift + (hz(j, zz) - 1) * hJprod(hnz - zz);
%     end
%     hz(j, hnz) = (cz_pr(j) - shift);
% end
% 
% % -------------------------------------------------------------------------
% 
% % Вычисление функции правдоподобия от истинных параметров иерахической
% % модели
% % hLreal = - HL(halpha, x, hz, hJ, N, hJprod, hJshift, hnz, halpha_dim)
% 
% % hLreal = - HL(halpha1, x, hz, hJ, N, hJprod, hJshift, hnz)
% 
% % Начальное приближение параметров и ЛФП от начального приближения
% % иерахической модели
% 
% % halpha0 = zeros(halpha_dim, Xdim);
% halpha0 = [0;0;0];
% % halpha0 = [0;0;0;0;0;0;0;0];
% % hL0 = - HL(halpha0, x, hz, hJ, N, hJprod, hJshift, hnz, halpha_dim)
% 
% % Поиск максимума функции правдоподобия иерахической модели
% i
% % options = optimset('LargeScale', 'off', 'GradObj', 'off', 'Display', 'iter', 'MaxIter', 5);
% % [halpha, hFVAL, hEXITFLAG, houtput, hgrad, hhessian] = fminunc('HL', halpha0, options, x, hz, hJ, N, hJprod, hJshift, hnz, halpha_dim);
% [halpha, hFVAL] = fminunc('HL', halpha0, options, x, hz, hJ, N, hJprod, hJshift, hnz, halpha_dim)
% quality_halpha = QualHZ(halpha, hreal, x, N, hJ, z, hJprod, hJshift, hnz)
% % halpha, hFVAL, hEXITFLAG, houtput, hgrad, hhessian
% 
% % % -------------------------------------------------------------------------
% % Перенумерация откликов
% temp = hz(:,1);
% for j = 1:hnz-1
%     hz(:, j) = hz(:, j + 1);
% end
% hz(:, hnz) = temp;
% clear temp;
% [halpha_wrong, hFVALwrong] = fminunc('HL', halpha0, options, x, hz, hJ, N, hJprod, hJshift, hnz, halpha_dim)
% quality_halpha_wrong = QualHZrenum(halpha_wrong, hreal, x, N, hJ, z, hJprod, hJshift, hnz)
% % % % -------------------------------------------------------------------------
% % 
% % % -------------------------------------------------------------------------
% % Перестроение векторов откликов для перекрестного оценивания
% cz_fr_hz = cz;
% hz_fr_cz = hz;
% for j = 1:N
%     shift = hz(j, hnz);
%     for zz = 1:hnz - 1
%         shift = shift + (hz(j, zz) - 1) * hJprod(hnz - zz);
%     end
%     cz_fr_hz(j) = shift;
%     shift = 0;
%     for zz = 1:hnz - 1
%        hz_fr_cz(j, zz) = ceil((cz(j) - shift) / hJprod(hnz - zz));
%        shift = shift + (hz_fr_cz(j, zz) - 1) * hJprod(hnz - zz);
%     end
%     hz_fr_cz(j, hnz) = (cz(j) - shift);
% end
% % 
% % % -------------------------------------------------------------------------
% % % Проверка правильности перекодирования откликов
% % % 
% % % cz_fr_hz2 = cz_fr_hz;
% % % hz_fr_cz2 = hz_fr_cz;
% % % for j = 1:N
% % %     shift = hz_fr_cz(j, hnz);
% % %     for zz = 1:hnz - 1
% % %         shift = shift + (hz_fr_cz(j, zz) - 1) * hJprod(hnz - zz);
% % %     end
% % %     cz_fr_hz2(j) = shift;
% % %     shift = 0;
% % %     for zz = 1:hnz - 1
% % %        hz_fr_cz2(j, zz) = ceil((cz_fr_hz(j) - shift) / hJprod(hnz - zz));
% % %        shift = shift + (hz_fr_cz2(j, zz) - 1) * hJprod(hnz - zz);
% % %     end
% % %     hz_fr_cz2(j, hnz) = (cz_fr_hz(j) - shift);
% % % end
% % 
% % % -------------------------------------------------------------------------
% % Оценивание параметров классической модели на откликах иерархической
% [calpha_onhz, cFVAL_onhz] = fminunc('CL', calpha0, options, x, cz_fr_hz, cJ, N)
% % -------------------------------------------------------------------------
% calpha_onhz = [calpha_onhz(1); calpha_onhz(3); calpha_onhz(2)]
% % -------------------------------------------------------------------------
% quality_calpha_onhz = QualCZ(calpha_onhz, creal, x, N, cJ, z(N, 1))
% % % -------------------------------------------------------------------------
% % % Оценивание параметров иерархической модели на откликах классической
% % % [halpha_oncz, hFVAL_oncz] = fminunc('HL', halpha0, options, x, hz_fr_cz, hJ, N, hJprod, hJshift, hnz, halpha_dim)
% % % quality_halpha_oncz = QualHZ(halpha_oncz, hreal, x, N, hJ, z, hJprod, hJshift, hnz)
% % % -------------------------------------------------------------------------
% % % Информационные критерии
% % 
% % % cAIC = 2 * cFVAL + 2 * size(calpha, 1)
% % % cBIC = 2 * cFVAL + log(N) * size(calpha, 1)
% % % cHQ = 2 * cFVAL + 2* log(log(N)) * size(calpha, 1)
% % % cKIC = 2 * cFVAL + 3 * size(calpha, 1)
% % % 
% % % hAIC = 2 * hFVAL + 2 * size(halpha, 1)
% % % hBIC = 2 * hFVAL + log(N) * size(halpha, 1)
% % % hHQ = 2 * hFVAL + 2* log(log(N)) * size(halpha, 1)
% % % hKIC = 2 * hFVAL + 3 * size(halpha, 1)
% % % 
% % % CRIT = [cAIC; cBIC; cHQ; cKIC; hAIC; hBIC; hHQ; hKIC];
% % % 
% % % [C, I] = min(CRIT)
% 
% -------------------------------------------------------------------------
% Накапливание оценок и поиск усредненных

alpha1(:, i) = calpha;
alphar(:, i) = ralpha;
alphaco(:, i) = coalpha;
% alpha2(:, i) = halpha;
% % alpha3(:, i) = calpha_onhz;
% % alpha4(:, i) = halpha_oncz;
% 
% Lval1(:, i) = cFVAL;
% Lval2(:, i) = hFVAL;
% % Lval3(:, i) = cFVAL_onhz;
% % Lval4(:, i) = hFVAL_oncz;
% 
% cest_shift(i) = sum((creal - calpha).^2);
% hest_shift(i) = sum((hreal - halpha).^2);
% % cest_shift_onhz(i) = sum((creal - calpha_onhz).^2);
% % hest_shift_oncz(i) = sum((hreal - halpha_oncz).^2);
% 
qual1(i) = quality_calpha;
% qualr(i) = quality_ralpha;
qualco(i) = quality_coalpha;
% qual2(i) = quality_halpha;
% % qual3(i) = quality_calpha_onhz;
% % qual4(i) = quality_halpha_oncz;

% -------------------------------------------------------------------------

% alpha1(:, i) = halpha;
% alpha2(:, i) = halpha_wrong;
% alpha3(:, i) = calpha_onhz;
% 
% Lval1(:, i) = hFVAL;
% Lval2(:, i) = hFVALwrong;
% Lval3(:, i) = cFVAL_onhz;
% 
% hest_shift(i) = sum((hreal - halpha).^2);
% hest_shift_wrong(i) = sum((hreal - halpha_wrong).^2);
% cest_shift_onhz(i) = sum((creal - calpha_onhz).^2);
% 
% qual1(i) = quality_halpha;
% qual2(i) = quality_halpha_wrong;
% qual3(i) = quality_calpha_onhz;


end

M1 = mean(alpha1, 2)
sum((creal - M1).^2)
Mr = mean(alphar, 2)
sum((rreal - Mr).^2)
Mco = mean(alphaco, 2)
sum((coreal - Mco).^2)
% sum((hreal - M1).^2)
% L1 = mean(Lval1, 2)
% M2 = mean(alpha2, 2)
% sum((hreal - M2).^2)
% L2 = mean(Lval2, 2)
% M3 = mean(alpha3, 2)
% sum((creal - M3).^2)
% L3 = mean(Lval3, 2)
% % M4 = mean(alpha4, 2)
% % sum((hreal - M4).^2)
% % L4 = mean(Lval4, 2)
% 
% % shift1 = mean(cest_shift)
% % shift2 = mean(hest_shift)
% 
% shift1 = mean(hest_shift)
% shift2 = mean(hest_shift_wrong)
% 
% % shift3 = mean(cest_shift_onhz)
% % shift4 = mean(hest_shift_oncz)
% 
mean(qual1)
% mean(qualr)
mean(qualco)
% mean(qual2)
% mean(qual3)
% % mean(qual4)

beep
% -------------------------------------------------------------------------