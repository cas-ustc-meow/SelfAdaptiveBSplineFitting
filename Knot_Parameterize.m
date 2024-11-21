function [knots,data] = Knot_Parameterize(datas,model,m,n)
%KNOT_PARAMETERIZE ������ڵ����²���������ɢ��Ϊ�ʺ���B�����Ľڵ��Ų�
%   ���룺datas��׼�������������ݵ���ż����Ϊx���꣬������Ϊy����
%   ���룺model��׼����ϵ�ģ�ͱ��
%   ���룺m��׼��ʹ�õ�B��������
%   ���룺n��׼��ʹ�õ�B��������������
%   �����knots������ʹ�õ�B�����ڵ����У���С��������
%   �����data��׼����ϵ����ݵ㣬���о��󣬵�һ����x���꣬�ڶ�����y���꣬�������ǻ���������t����

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

