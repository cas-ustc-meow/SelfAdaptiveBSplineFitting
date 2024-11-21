function [error,cx,cy,knots] = SelfAdaptedSplines_fit(m,n,datas,model,tolerance)
%SELFADAPTEDSPLINES_FIT 使用自适应B样条（节点插入算法）实现数据点拟合任务
%   输入：m：使用B样条的次数
%   输入：n：初始样条基函数个数
%   输入：datas：使用的拟合数据点表格，偶数行为x坐标，奇数行为y坐标
%   输入：model：准备测试的模型
%   输入：tolerance：允许的拟合误差，这里是最大模范数

%% 第一步：先行存储需要的节点和数据点信息
[knots,data]=Knot_Parameterize(datas,model,m,n);
num_knots=size(knots,1);
num_bases=num_knots-m-1;
num_data=size(data,1);
%% 第二步：对数据点进行反复拟合和加密，直到误差达到容许阈值
A=zeros(num_data,num_bases);%B样条拟合使用的法方程
x=data(:,1);
y=data(:,2);%B样条拟合使用的数据列向量
t=data(:,3);
error=1;
while error>=tolerance
    %使用deBoor-Cox算法计算拟合矩阵
    for i=1:num_data
        for j=1:num_bases
            A(i,j)=deBoor_Cox_bases(knots,j,data(i,3),m);
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
    [error,unacceptable]=fit_evaluation(cx,cy,knots,data,m,tolerance);
    if error_old<=error
        break;
    end
    %如果不在容许阈值以内，查找不在阈值以内的参数所在区间，并进行B样条基函数加密
    if error>=tolerance
        unacceptable(:)=1;
        unacceptable(1:m)=0;
        unacceptable(num_bases+1:m+num_bases)=0;
        new_knots=subdivide(knots,unacceptable);
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
    xsamples(1,i)=deBoor_Cox_algorithm(knots,cx,tsamples(1,i),m);
    ysamples(1,i)=deBoor_Cox_algorithm(knots,cy,tsamples(1,i),m);
end
plot(xsamples,ysamples,'LineWidth',3);
end

