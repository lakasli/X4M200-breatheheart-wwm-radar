function imf = emd(x)
% Empiricial Mode Decomposition (Hilbert-Huang Transform)
% EMD分解或HHT变换
% 返回值为cell类型，依次为一次IMF、二次IMF、...、最后残差
x   = transpose(x(:));
imf = [];
while ~ismonotonic(x)
    x1 = x;
    sd = Inf;
    while (sd > 0.1) || ~isimf(x1)
        s1 = getspline(x1);         % 极大值点样条曲线
        s2 = -getspline(-x1);       % 极小值点样条曲线
        x2 = x1-(s1+s2)/2;
       
        sd = sum((x1-x2).^2)/sum(x1.^2);
        x1 = x2;
    end
   
    imf{end+1} = x1;
    x          = x-x1;
end
imf{end+1} = x;
% 是否单调
function u = ismonotonic(x)
u1 = length(findpeaks(x))*length(findpeaks(-x));
if u1 > 0
    u = 0;
else
    u = 1;
end
% 是否IMF分量
function u = isimf(x)
N  = length(x);
u1 = sum(x(1:N-1).*x(2:N) < 0);                     % 过零点的个数
u2 = length(findpeaks(x))+length(findpeaks(-x));    % 极值点的个数
if abs(u1-u2) > 1
    u = 0;
else
    u = 1;
end
% 根据极大值点构造样条曲线
function s = getspline(x)
N = length(x);
[p, locs] = findpeaks(x);
s = spline([0, locs, N+1], [0, x(locs), 0], 1:N);