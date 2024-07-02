function [x_position,y_position,x_velocity,y_velocity]= initialize(L,N)
positions=zeros(N,2);%Empty matrix for positions
particle=0;%Counter for particle number
%Giving positions as cubic lattice points
for i=1:nthroot(N,2)
    for j=1:nthroot(N,2)
        particle=particle+1;
        positions(particle,:)=[(L*i)/nthroot(N,2)-1,(L*j)/nthroot(N,2)-1];
    end
end
%Extracting x&y values
x_position=positions(:,1);
y_position=positions(:,2);

%Giving random values of velocity between -1 & 1
x_velocity=(2*rand(N,1)-1);
y_velocity=(2*rand(N,1)-1);
%Shifting the velocity to make average 0
vx_avg=sum(x_velocity)/N;
vy_avg=sum(y_velocity)/N;
for i = 1:N
    x_velocity(i)=x_velocity(i)-vx_avg;
    y_velocity(i)=y_velocity(i)-vy_avg;
end

