% ���������� ���

%   ������� ���������:
% alpha - ������ ����������� ����������
% x - ����������
% z - ������ �������
% J - ����������� ������� �������
% N - ���������� �������������

%   �������� ���������:
% out - ������������� �������� ������� �������������

% -------------------------------------------------------------------------

% function [out, gradL] = HL(alpha0, x, hz, hJ, N, hnz)
function [out, gradL] = HL(alpha, x, hz, hJ, N, hJprod, hJshift, hnz, halpha_dim)
out2 = 0;
grad2 = zeros(halpha_dim, 1);
grad1 = grad2;

% ���������� ��� � ���������

for i = 1:N
    out1 = 0;

    % ���������� ������������ ��� ������� �������
    for j = 1:hJ(1)
%             ksij = zeros(halpha_dim, 1);
            P(1, j) = 1;
            subsum = 1;
%             subsumksi = 0;
            for k = 1:hJ(1)-1
%                 ksis = zeros(halpha_dim, 1);
                subsum = subsum + exp(x(i, :) * alpha(k, :));
%                 subsumksi = subsum + exp(x(i, :) * alpha(k, :)) * ksis;
            end
            if j == hJ(1)
                pr = 1;
            else
                pr =  exp(x(i, :) * alpha(j, :));
%                 ksij(j) = x(i, :);
%                 grad1 = grad1 + (ksij + ksij * (subsum - 1) - subsumksi) / (subsum);
            end
            P(1, j) = P(1, j) * pr / subsum;
            
    end

% A1 = exp(x(i, :) * alpha(1)) / (1 + exp(x(i, :) * alpha(1)));
% A2 = exp(x(i, :) * alpha(1)) / (1 + exp(x(i, :) * alpha(1)));
% A3 = 1 / (1 + exp(x(i, :) * alpha(1)));
% A4 = 1 / (1 + exp(x(i, :) * alpha(1)));
% 
% P(1,1)=A1;
% P(1,2)=A2;
% P(1,3)=A3;
% P(1,4)=A4;

    out1 = out1 + log(P(1, hz(i, 1)));

% -------------------------------------------------------------------------
    for k = 2:hnz                               % ���� �� ��������
        Ptmp = ones(1, hJ(k));                  % ����������� �������
        outshift = hz(i, k - 1);                % ����� �� ����������
        for kk = 3:k
            outshift = outshift + (hz(i, kk - 1) - 1) * hJprod(kk - 2);
        end
%         grad2 = 0;
        for j = 1:hJ(k)                         % ���� �� ���������
%             ksij = zeros(halpha_dim, 1);
            % ���������� ������������
            subsum = 1;
%             subsumksi = 0;
            for s = 1:hJ(k)-1
%                 ksis = zeros(halpha_dim, 1);
%                 ksis(hJshift(k) + outshift + hJprod(k - 1) * (s - 1)) = x(i, :);
                subsum = subsum + exp(x(i, :) * alpha(hJshift(k) + outshift + hJprod(k - 1) * (s - 1), :));
%                 subsumksi = subsumksi + exp(x(i, :) * alpha(hJshift(k) + outshift + hJprod(k - 1) * (s - 1), :)) * ksis;
            end
            if j == hJ(k)
                pr = 1;
            else
                pr =  exp(x(i, :) * alpha(hJshift(k) + outshift + hJprod(k - 1) * (j - 1), :));
                
%                 ksij(hJshift(k) + outshift + hJprod(k - 1) * (j - 1)) = x(i, :);
%                 grad2 = grad2 + (ksij + ksij * (subsum - 1) - subsumksi) / (subsum);
            end
            Ptmp(j) = Ptmp(j) * pr / subsum;
                        
        end
        P(k, :) = Ptmp;

% A1 = exp(x(i, :) * alpha(2)) / (1 + exp(x(i, :) * alpha(2)));
% A2 = 1 / (1 + exp(x(i, :) * alpha(2)));
% A3 = exp(x(i, :) * alpha(3)) / (1 + exp(x(i, :) * alpha(3)));
% A4 = 1 / (1 + exp(x(i, :) * alpha(3)));
% 
% P(2,1)=A1;
% P(2,2)=A2;
% P(2,3)=A3;
% P(2,4)=A4;

        out1 = out1 + log(P(k, hz(i, k)));
%         grad1 = grad1 + grad2;
    end
    out2 = out2 + out1;
end

out = - out2;
gradL = - grad1;