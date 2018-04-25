function [position,elementPointX,elementPointY] = array_calculate_bak(a,b,num_of_element)  
%参数意义
%a,长轴半径
%b,短轴半径
% num_of_element  椭圆上要放置的阵元个数
%position,椭圆阵列中阵元的极坐标位置
%k 长轴为a，短轴为b的椭圆上阵元的个数为 k
%elementPointX,阵元在直角坐标系下的X轴的坐标位置
%elementPointY，阵元在直角坐标系下的Y轴的坐标位置
% a=4,b=3,num_of_element=10;
minSpace=(2*pi*b+4*(a-b))/num_of_element;
h=0.000001; %步长
data=0:h:a;  
x1=fliplr(data);  %将data进行反转
y1=sqrt(b^2-b^2*x1.^2/a^2); %椭圆公式
x2=-data;
y2=sqrt(b^2-b^2*x2.^2/a^2);
x3=-fliplr(data);
y3=-sqrt(b^2-b^2*x3.^2/a^2);
x4=data;
y4=-sqrt(b^2-b^2*x4.^2/a^2);
%y4=-fliplr(y1);
x=[x1,x2,x3,x4];
y=[y1,y2,y3,y4];
n=length(x);  %份数
lamda=minSpace;   %lamda波长
i=1;  %第i个点
k=1;   %没有算初始点，最终导致k=num_of_element+1,
% e=0.0005;   %误差,要是步长的n倍，否则会导致某个点跳过之后，就一直没有合适的点
sumlength=0;
elementPointX=[];  
elementPointY=[];
%while (sum<2*pi*b+4*(a-b))
elementPointX(k)=a;
elementPointY(k)=0;
while(1<2)
    while (i<n)  %n点的个数
        m=sqrt( (x(i+1)-x(i))^2+(y(i+1)-y(i))^2);  %两点之间的距离,m突然一下变的很大，不在lamda范围内了
        sumlength=sumlength+m;
        if((lamda-lamda*0.00001)<sumlength && sumlength<=(lamda+lamda*0.00001))  %点i在误差范围内时，保存
               k=k+1;
    %            v(k)=i;        
               sumlength=0;
               elementPointX(k)=x(i+1);
               elementPointY(k)=y(i+1);
               i = i + 1;
         else
                i = i + 1;
         end
    end 
    if(lamda-lamda*0.05<sumlength)
        break;
    else
        lamda=lamda-lamda*0.001;
        i=1;  
        k=1;   
        sumlength=0;
    end
end
%plot(elementPointX,elementPointY);  //直角坐标

angle=0;%角度
radius=0;%半径
len=int8(length(elementPointX));
position=zeros(len,2); %存放椭圆上面的点的半径和角度

for index=1:len
    radius=sqrt(elementPointX(index)^2+elementPointY(index)^2);  %半径
    if (elementPointX(index)>0 && elementPointY(index)>0)
        angle=atan(elementPointY(index)/elementPointX(index));
    else if (elementPointX(index)<0 && elementPointY(index)>0)
        angle=pi/2+abs(atan(elementPointX(index)/elementPointY(index)));
    else if (elementPointX(index)<0 && elementPointY(index)<0)
        angle=pi+abs(atan(elementPointY(index)/elementPointX(index)));
    else if (elementPointX(index)>0 && elementPointY(index)<0) 
        angle=pi*3/2+abs(atan(elementPointX(index)/elementPointY(index)));
        end
        end
        end
    end
    position(index,1)=radius;  %半径
    position(index,2)=angle;   %角度
end
% position


