function [knots,data] = Knot_Parameterize(datas,model,m,n)
%KNOT_PARAMETERIZE 对输入节点重新参数化和离散成为适合于B样条的节点排布
%   输入：datas：准备导入的拟合数据点表格，偶数行为x坐标，奇数行为y坐标
%   输入：model：准备拟合的模型编号
%   输入：m：准备使用的B样条次数
%   输入：n：准备使用的B样条基函数个数
%   输出：knots：初步使用的B样条节点序列，从小到大排序
%   输出：data：准备拟合的数据点，三列矩阵，第一列是x坐标，第二列是y坐标，第三列是弧长参数化t坐标

data=csvread(datas);
data=data(2*model:2*model+1,:)';
%data=sortrows(data',1);

t_coords=zeros(size(data,1),1);
for i=2:size(data,1)
    t_coords(i,1)=t_coords(i-1)+sqrt((data(i,1)-data(i-1,1))^2+(data(i,2)-data(i-1,2))^2);
end
t_coords=t_coords/(t_coords(end,1));
data=[data t_coords];

knots=zeros(n+m+1,1);
x_min=min(data(:,3));
x_max=max(data(:,3));
knot_max=x_max;
knot_min=x_min;
knots(1:m+1)=knot_min;
knots(n+1:n+m+1)=knot_max;
temp=linspace(knot_min,knot_max,n-m+1);
knots(m+2:n)=temp(2:n-m);
% knots=linspace(knot_min,knot_max,16+m+1)';
end

