function [dx,dy,r]=separation(N,L,x_position,y_position)
dx=zeros(N,1);
dy=zeros(N,1);

r=zeros(N,1);
for i=1:N
    for j= 1:N
    dx(i,j)=x_position(i)-x_position(j);
    if abs(dx(i,j))>L/2
        dx(i,j)=(L-abs(dx(i,j)))*(-dx(i,j))/abs(dx(i,j));
    end
    dy(i,j)=y_position(i)-y_position(j);
    if abs(dy(i,j))>L/2
        dy(i,j)=(L-abs(dy(i,j)))*(-dy(i,j))/abs(dy(i,j));
    end
    
    r(i,j)=sqrt(dx(i,j)^2+dy(i,j)^2);
    end
end
