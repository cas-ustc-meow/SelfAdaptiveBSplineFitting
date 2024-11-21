function c = NormalEquation(A,y)
%NORMALEQUATION 求解拟合需要的法方程
%   输入：A：法矩阵
%   输入：y：拟合数值
%   输出：c：拟合系数
c=(A'*A)\(A'*y);
end

