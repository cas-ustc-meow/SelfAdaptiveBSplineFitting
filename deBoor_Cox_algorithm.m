function result = deBoor_Cox_algorithm(knots,c,sample,m)
%DEBOOR_COX_ALGORITHM 使用deBoor-Cox算法实现B样条的函数值计算，用于绘图和估计
%   输入：knots：使用B样条基函数的节点
%   输入：c：使用B样条函数的组合系数
%   输入：sample：准备估计的采样点
%   输入：m：使用B样条的次数
%   输出：result：B样条拟合函数在当前位置的数值

%先找出当前采样点所在的节点区间编号
niu=0;
for i=1:size(knots,1)-1
    if sample>=knots(i) && sample<knots(i+1)
        niu=i;
        break;
    else
        niu=size(knots,1)-m-1;
    end
end
%利用deBoor-Cox迭代计算基函数取值——此时只有对应基函数权重为1，其余全部为0
min_base=niu-m;
max_base=niu;
c_tower=zeros(m+1,m+1);
for i=1:m+1
    c_tower(1,i)=c(min_base+i-1);
end
for r=1:m
    for j=r+1:m+1
        knot_number=min_base+j-1;
        c_tower(r+1,j)=(knots(knot_number+m+1-r)-sample)/(knots(knot_number+m+1-r)-knots(knot_number))*c_tower(r,j-1)+(sample-knots(knot_number))/(knots(knot_number+m+1-r)-knots(knot_number))*c_tower(r,j);
    end
end
result=c_tower(m+1,m+1);

end

