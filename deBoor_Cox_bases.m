function result = deBoor_Cox_bases(knots,base_number,sample,m)
%DEBOOR_COX_ALGORITHM 使用deBoor-Cox算法实现B样条基函数的函数值计算，用于构建法方程
%   输入：knots：使用B样条基函数的节点序列
%   输入：base_number：准备参加估计的基函数编号
%   输入：sample：准备估计的采样点
%   输入：m：使用B样条的次数
%   输出：result：B样条基函数在当前位置的数值

%% 判断采样点是否在当前B样条基函数的紧支集下，如果不是就自动返回0；如果是则继续计算
if sample<=knots(base_number) || sample>=knots(base_number+m+1)
    result=0;
else
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
        if min_base+i-1==base_number
            c_tower(1,i)=1;
        else
            c_tower(1,i)=0;
        end
    end
    for r=1:m
        for j=r+1:m+1
            knot_number=min_base+j-1;
            c_tower(r+1,j)=(knots(knot_number+m+1-r)-sample)/(knots(knot_number+m+1-r)-knots(knot_number))*c_tower(r,j-1)+(sample-knots(knot_number))/(knots(knot_number+m+1-r)-knots(knot_number))*c_tower(r,j);
        end
    end
    result=c_tower(m+1,m+1);
end
end

