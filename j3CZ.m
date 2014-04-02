% ���������� ������� ��������

%     ������� ���������:
% alpha - ������ ����������� ����������
% x - ����������
% N - ���������� �������������
% J - ���������� �������� �������
% z - ��������� ������ ��������

%     �������� ���������:
% z - ��������������� ������ ��������

% -------------------------------------------------------------------------

function out = j3CZ(alpha, xx, N, J, z, lambda, eps)

P = ones(3, 1);
out = zeros(N, 1);
flag = 0;

for i = 1:N
    x = xx(i,1);
    % ���������� ������������ ��������� � ��������
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

% P(1) = exp(x(i, :) * alpha(1, :)) / (1 + exp(x(i, :) * alpha(1, :)) + exp(x(i, :) * alpha(2, :)) + exp(x(i, :) * alpha(3, :)));
% P(2) = exp(x(i, :) * alpha(2, :)) / (1 + exp(x(i, :) * alpha(1, :)) + exp(x(i, :) * alpha(2, :)) + exp(x(i, :) * alpha(3, :)));
% P(3) = exp(x(i, :) * alpha(3, :)) / (1 + exp(x(i, :) * alpha(1, :)) + exp(x(i, :) * alpha(2, :)) + exp(x(i, :) * alpha(3, :)));
% P(4) = 1 / (1 + exp(x(i, :) * alpha(1, :)) + exp(x(i, :) * alpha(2, :)) +
% exp(x(i, :) * alpha(3, :)));

P(1) = exp(alpha(1,1) + x * alpha(2,1) + x * x * alpha(3,1)) / (1 + exp(alpha(1,1) + x * alpha(2,1) + x * x * alpha(3,1)) + exp(alpha(1,2) + x * alpha(2,2) + x * x * alpha(3,2)));
P(2) = exp(alpha(1,2) + x * alpha(2,2) + x * x * alpha(3,2)) / (1 + exp(alpha(1,1) + x * alpha(2,1) + x * x * alpha(3,1)) + exp(alpha(1,2) + x * alpha(2,2) + x * x * alpha(3,2)));
P(3) = 1 - P(1) - P(2);

% P(1) = exp(alpha(1, 1) + x(i, :) * alpha(2, 1) + x(i, :) * x(i, :) * alpha(3, 1))...
%     / (1 + exp(alpha(1, 1) + x(i, :) * alpha(2, 1) + x(i, :) * x(i, :) * alpha(3, 1))...
%     + exp(alpha(1, 2) + x(i, :) * alpha(2, 2) + x(i, :) * x(i, :) * alpha(3, 2)));
% 
% P(2) = exp(alpha(1, 2) + x(i, :) * alpha(2, 2) + (x(i, :)^2) * alpha(3, 2))...
%     / (1 + exp(alpha(1, 1) + x(i, :) * alpha(2, 1) + (x(i, :)^2) * alpha(3, 1))...
%     + exp(alpha(1, 2) + x(i, :) * alpha(2, 2) + (x(i, :)^2) * alpha(3, 2)));
% 
% P(3) = 1 - P(1) - P(2)
% i
% x(i,:)

delta = 0;
for m = 1:J-1
    delta = delta + (P(m) ^ (1 - lambda));
end

cP(1) = (P(1) ^ (1 - lambda)) / delta;
cP(2) = (P(2) ^ (1 - lambda)) / delta;
cP(3) = (P(3) ^ (1 - lambda)) / delta;

cP(1) = (1 - P(1)) / (J - 1);
cP(2) = (1 - P(2)) / (J - 1);
cP(3) = (1 - P(3)) / (J - 1);

% ���������

% for m=1:J
%     P(m) = P(m) * (1 - eps) + eps * cP(m);
% end

    % ���������� ������������ ������������
    for j = 1:2
        P(j + 1) = P(j) + P(j + 1);
    end
    
    if P(3) ~= 1
        P(3) = 1;
    end

    % ���������� �������� ������� � ����������� �� ��������� � �����������
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
    out(i)
    flag = 0;
end

end