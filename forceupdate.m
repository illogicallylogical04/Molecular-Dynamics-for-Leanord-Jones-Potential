function [Total_Force_x,Total_Force_y]=forceupdate(N,Fc,dx,dy,r)
%Creating empty matrices for the components of force
rc=5;
Forces_x=zeros(N,N);
Forces_y=zeros(N,N);
for i=1:N
    for j= 1:N
        %Implementing the cut-off distance
        if r(i,j)==0 || r(i,j)>rc
           Forces_x(i,j)=0;Forces_y(i,j)=0;
        else
           %Calculating the components of force for each pair of particles
           Forces_x(i,j)=(96*((1/(r(i,j)))^13-((1/2)*(1/(r(i,j))^7)))*(dx(i,j)/r(i,j)))-Fc;
           Forces_y(i,j)=(96*((1/(r(i,j)))^13-((1/2)*(1/(r(i,j))^7)))*(dy(i,j)/r(i,j)))-Fc;
        end
    end
end
%Summing up the rows to find out the total force on each individual particle
Total_Force_x=(sum(Forces_x))';
Total_Force_y=(sum(Forces_y))';
end


