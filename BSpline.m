%作业6：实现简单的三次B样条
figure; 

%%
%首先完成取点
h = drawpolyline;
hold on;
p=h.Position;
n=size(p,1);n=n-1;

%%
%其次先行准备参数节点，由于规定4阶B样条，故有n+4个节点
%为了方便计算，目前使用整数节点
global k;
k=4;
global knots;
knots=zeros(n+k+1,1);
for i=k+1:n+2
    knots(i,1)=i-k;
end
knots(1)=0;knots(2)=0;knots(3)=0;
knots(n+k+1)=n+2-k;knots(n+k)=n+2-k;knots(n+k-1)=n+2-k;
%%
%构造B样条de Boor算法需要的塔式结构
global K;
K=zeros(k,k,2);

%%
global t;
t=(0:0.001*(n+2-k):n+2-k);
%re=calculateBezier(h.Position, t);
%hcurve = plot(calculateBezier(h.Position, t), 'g', 'linewidth', 2);
%hcurve=calculateBezier(h.Position, t,hcurve);
re=calculateBSpline(h.Position);
hcurve = plot(re(:,1),re(:,2), 'g', 'linewidth', 2);
h.addlistener('MovingROI', @(h, evt) calculateBSpline(evt.CurrentPosition, hcurve));

%%
%递归形式的de Boor算法
function re=calculateBSpline(p,h)
global K;
global t;
global k;
global knots;
    re=[];
    n=size(p,1)-1;
    
    for i=1:1000
        %先查找到当前t所在的区间
        r=floor(t(i))+k-1;
        %装填de Boor的迭代初值
        for s=1:k
            K(s,1,1)=p(r-k+1+s,1);
            K(s,1,2)=p(r-k+1+s,2);
        end
        %迭代出所求参数处的B样条取值
        for j=1:k-1
            for m=r-k+1+j:r
                if knots(m+k-j+1)-knots(m+1)==0
                    alpha=0;
                else
                    alpha=(t(i)-knots(m+1))/(knots(m+k-j+1)-knots(m+1));
                end
                K(m-r+k,j+1,1)=(1-alpha)*K(m-r+k-1,j,1)+alpha*K(m-r+k,j,1);
                K(m-r+k,j+1,2)=(1-alpha)*K(m-r+k-1,j,2)+alpha*K(m-r+k,j,2);
            end
        end
        %取出需要的结果
        re=[re;K(k,k,1),K(k,k,2)];
    end
    re=[re;p(n+1,1),p(n+1,2)];
    % evaluate the function at x2
    if nargin>1,set(h,'xdata',re(:,1),'ydata',re(:,2));%关键问题在这里！
    end
end