function [error,cx,cy,knots] = HBSplines_fit(m,n,datas,model,tolerance,hierarchy)
%HBSPLINES_FIT ʹ�ò��B����ʵ�����ݵ��������
%   ���룺m��ʹ��B�����Ĵ���
%   ���룺n����ʼ��������������
%   ���룺datas��ʹ�õ�������ݵ���ż����Ϊx���꣬������Ϊy����
%   ���룺model��׼�����Ե�ģ��
%   ���룺tolerance���������������������ģ����
%   ���룺hierarchy�����B�������趨����

%% ��һ�������д洢��Ҫ�Ľڵ�����ݵ���Ϣ
[knots,data]=Knot_Parameterize(datas,model,m,n);
num_knots=size(knots,1);
num_bases=num_knots-m-1;
num_data=size(data,1);
hierarchy_basis=zeros(num_bases,1);
hierarchy_knots=zeros(num_knots,1);
%% �ڶ����������ݵ���з�����Ϻͼ��ܣ�ֱ�����ﵽ������ֵ���߲�δﵽָ������
A=zeros(num_data,num_bases);%B�������ʹ�õķ�����
x=data(:,1);
y=data(:,2);%B�������ʹ�õ�����������
t=data(:,3);
error=1;
current_level=0;
while error>=tolerance && current_level<=hierarchy
    %ʹ��deBoor-Cox�㷨������Ͼ���ע��Ҫ������B����
    for i=1:num_data
        for j=1:num_bases
            A(i,j)=deBoor_Cox_bases_HB(knots,j,data(i,3),m,hierarchy_basis,hierarchy_knots);
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
    [error,unacceptable]=fit_evaluation_HB(cx,cy,knots,data,m,tolerance,hierarchy_basis,hierarchy_knots);
    if error_old<=error
        break;
    end
    %�������������ֵ���ڣ����Ҳ�����ֵ���ڵĲ����������䣬�����в��B����������ϸ��
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

%% �����������ӻ����յ���Ͻ����Ŀ�����ݵ�
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

