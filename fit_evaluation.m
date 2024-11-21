function [error,unacceptable] = fit_evaluation(cx,cy,knots,data,m,tolerance)
%FIT_EVALUATION 评估B样条拟合的效果是否在给定阈值以内
%   输入：c：当前B样条函数的组合系数
%   输入：knots：当前B样条基函数的节点序列
%   输入：data：目标拟合的数据点x-y坐标矩阵
%   输入：m：使用的B样条次数
%   输入：tolerance：允许的拟合误差，这里是最大模范数
%   输出：error：B样条拟合给定数据点的误差，目前使用最大模范数
%   输出：unacceptable：当前B样条节点序列中存在误差超过阈值的区间指示，1表示有违反，0表示无违反

error_list=zeros(size(data,1),1);
for i=1:size(data,1)
    fit_value_x=deBoor_Cox_algorithm(knots,cx,data(i,3),m);
    fit_value_y=deBoor_Cox_algorithm(knots,cy,data(i,3),m);
    error_list(i,1)=sqrt((fit_value_x-data(i,1))^2+(fit_value_y-data(i,2))^2);
end

error=max(error_list);
unacceptable=zeros(size(knots,1)-1,1);
for j=1:size(data,1)
    if error_list(j,1)>=tolerance
        for i=1:size(knots,1)-1
            if data(j,1)>=knots(i) && data(j,1)<knots(i+1)
                unacceptable(i)=1;
                break;
            end
        end
    end
end


end

