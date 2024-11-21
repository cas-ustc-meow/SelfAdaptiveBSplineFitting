function [error,cx,cy,knots] = SelfAdaptedSplines_fit(m,n,datas,model,tolerance)
%SELFADAPTEDSPLINES_FIT ʹ������ӦB�������ڵ�����㷨��ʵ�����ݵ��������
%   ���룺m��ʹ��B�����Ĵ���
%   ���룺n����ʼ��������������
%   ���룺datas��ʹ�õ�������ݵ���ż����Ϊx���꣬������Ϊy����
%   ���룺model��׼�����Ե�ģ��
%   ���룺tolerance���������������������ģ����

%% ��һ�������д洢��Ҫ�Ľڵ�����ݵ���Ϣ
[knots,data]=Knot_Parameterize(datas,model,m,n);
num_knots=size(knots,1);
num_bases=num_knots-m-1;
num_data=size(data,1);
%% �ڶ����������ݵ���з�����Ϻͼ��ܣ�ֱ�����ﵽ������ֵ
A=zeros(num_data,num_bases);%B�������ʹ�õķ�����
x=data(:,1);
y=data(:,2);%B�������ʹ�õ�����������
t=data(:,3);
error=1;
while error>=tolerance
    %ʹ��deBoor-Cox�㷨������Ͼ���
    for i=1:num_data
        for j=1:num_bases
            A(i,j)=deBoor_Cox_bases(knots,j,data(i,3),m);
        end
    end
    %�����Ӧ�����̵õ�����õ�B����ϵ��
%     if rank(A)<num_bases
%         break;
%     end
    cx=NormalEquation(A,x);
    cy=NormalEquation(A,y);
    %�����������Ƿ���������ֵ����
    error_old=error;
    [error,unacceptable]=fit_evaluation(cx,cy,knots,data,m,tolerance);
    if error_old<=error
        break;
    end
    %�������������ֵ���ڣ����Ҳ�����ֵ���ڵĲ����������䣬������B��������������
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

%% �����������ӻ����յ���Ͻ����Ŀ�����ݵ�
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

