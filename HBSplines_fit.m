function [error,cx,cy,knots] = HBSplines_fit(m,n,datas,model,tolerance,hierarchy)
%HBSPLINES_FIT 使用层次B样条实现数据点拟合任务
%   输入：m：使用B样条的次数
%   输入：n：初始样条基函数个数
%   输入：datas：使用的拟合数据点表格，偶数行为x坐标，奇数行为y坐标
%   输入：model：准备测试的模型
%   输入：tolerance：允许的拟合误差，这里是最大模范数
%   输入：hierarchy：层次B样条的设定层数

%% 第一步：先行存储需要的节点和数据点信息
[knots,data]=Knot_Parameterize(datas,model,m,n);
num_knots=size(knots,1);
num_bases=num_knots-m-1;
num_data=size(data,1);
hierarchy_basis=zeros(num_bases,1);
hierarchy_knots=zeros(num_knots,1);
%% 第二步：对数据点进行反复拟合和加密，直到误差达到容许阈值或者层次达到指定上限
A=zeros(num_data,num_bases);%B样条拟合使用的法方程
x=data(:,1);
y=data(:,2);%B样条拟合使用的数据列向量
t=data(:,3);
error=1;
current_level=0;
while error>=tolerance && current_level<=hierarchy
    %使用deBoor-Cox算法计算拟合矩阵，注意要加入层次B样条
    for i=1:num_data
        for j=1:num_bases
            A(i,j)=deBoor_Cox_bases_HB(knots,j,data(i,3),m,hierarchy_basis,hierarchy_knots);
        end
    end
    %计算对应法方程得到拟合用的B样条系数
%     if rank(A)<num_bases
%         break;
%     end
    cx=NormalEquation(A,x);
    cy=NormalEquation(A,y);
    %计算拟合误差是否在容许阈值以内
    error_old=error;
    [error,unacceptable]=fit_evaluation_HB(cx,cy,knots,data,m,tolerance,hierarchy_basis,hierarchy_knots);
    if error_old<=error
        break;
    end
    %如果不在容许阈值以内，查找不在阈值以内的参数所在区间，并进行层次B样条基函数细分
    if error>=tolerance
        unacceptable(:)=1;
        unacceptable(1:m)=0;
        unacceptable(num_bases+1:m+num_bases)=0;
        new_knots=HB_subdivide(knots,unacceptable,m);
        if size(new_knots,1)-m-1<num_data
            knots=new_knots;
            num_knots=size(knots,1);
            num_bases=num_knots-m-1;
            A=zeros(num_data,num_bases);
        else
            break;
        end
    end
end

%% 第三步：可视化最终的拟合结果和目标数据点
figure(1);
scatter(data(:,1),data(:,2));
hold on;
tsamples=linspace(min(data(:,3)),max(data(:,3)),1000);
xsamples=zeros(1,size(tsamples,2));
ysamples=zeros(1,size(tsamples,2));
for i=1:size(tsamples,2)
    xsamples(1,i)=deBoor_Cox_algorithm_HB(knots,cx,tsamples(1,i),m);
    ysamples(1,i)=deBoor_Cox_algorithm_HB(knots,cy,tsamples(1,i),m);
end
plot(xsamples,ysamples,'LineWidth',3);
end

