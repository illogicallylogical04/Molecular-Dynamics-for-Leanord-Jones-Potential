function [Kinetic_Energy,Potential_Energy,Total_Energy]=Energycalc(N,v,r,rc,Fc,Vc)
%Total Kinetic energy of the system
Kinetic_Energy=(sum(v.^2))/(2*N);

%Creating empty matrix for potential energy values
Potential=zeros(N,N);
%Calculating the potential energy by considering each pair only once
for i=1:N-1
    for j= i+1:N
        if r(i,j)<=rc
           Potential(i,j)=8*((1/r(i,j)^12)-(1/r(i,j)^6))+Fc*r(i,j)-Fc*rc-Vc;
        end
    end
end
Potential_Energy=abs(sum(sum(Potential)))/N;%Total Potential Energy
Total_Energy=Kinetic_Energy+Potential_Energy;%Total Energy

