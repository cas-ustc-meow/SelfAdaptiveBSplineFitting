function result = deBoor_Cox_bases(knots,base_number,sample,m)
%DEBOOR_COX_ALGORITHM ʹ��deBoor-Cox�㷨ʵ��B�����������ĺ���ֵ���㣬���ڹ���������
%   ���룺knots��ʹ��B�����������Ľڵ�����
%   ���룺base_number��׼���μӹ��ƵĻ��������
%   ���룺sample��׼�����ƵĲ�����
%   ���룺m��ʹ��B�����Ĵ���
%   �����result��B�����������ڵ�ǰλ�õ���ֵ

%% �жϲ������Ƿ��ڵ�ǰB�����������Ľ�֧���£�������Ǿ��Զ�����0����������������
if sample<=knots(base_number) || sample>=knots(base_number+m+1)
    result=0;
else
    %���ҳ���ǰ���������ڵĽڵ�������
    niu=0;
    for i=1:size(knots,1)-1
        if sample>=knots(i) && sample<knots(i+1)
            niu=i;
            break;
        else
            niu=size(knots,1)-m-1;
        end
    end
    %����deBoor-Cox�������������ȡֵ������ʱֻ�ж�Ӧ������Ȩ��Ϊ1������ȫ��Ϊ0
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

