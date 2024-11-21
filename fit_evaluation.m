function [error,unacceptable] = fit_evaluation(cx,cy,knots,data,m,tolerance)
%FIT_EVALUATION ����B������ϵ�Ч���Ƿ��ڸ�����ֵ����
%   ���룺c����ǰB�������������ϵ��
%   ���룺knots����ǰB�����������Ľڵ�����
%   ���룺data��Ŀ����ϵ����ݵ�x-y�������
%   ���룺m��ʹ�õ�B��������
%   ���룺tolerance���������������������ģ����
%   �����error��B������ϸ������ݵ����Ŀǰʹ�����ģ����
%   �����unacceptable����ǰB�����ڵ������д���������ֵ������ָʾ��1��ʾ��Υ����0��ʾ��Υ��

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

