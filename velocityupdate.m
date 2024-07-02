function[v_new]=velocityupdate(N,velocity,dt,Fo,Fn)
%Create empty matrix for velocity of each particle
v_new=zeros(1,N);
for i=1:N
    v_new(i)=velocity(i)+(dt/2)*(Fo(i)+Fn(i));
end
%v_new(i) is the updated velocity of the ith particle
%velocity(i) is the velocity at  previous timestep of the ith particle
%dt is delta t
%Fo(i) is the force at the previous timestep of the ith particle
%Fn(i) is the force at the current timestep of the ith particle

