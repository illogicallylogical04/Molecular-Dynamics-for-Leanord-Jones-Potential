clc
clear
%% Defining the variables
dt=0.005;
dt2=dt*dt;
timesteps=2000;
tmax=dt*timesteps;
t=(0:dt:tmax);
N=100;
L=5*sqrt(N);
rc=5;
Fc=96*((1/(rc))^13-((1/2)*(1/(rc)^7)));
Vc=8*((1/rc)^12-(1/rc)^6);
%% initializing
[x_position,y_position,x_velocity,y_velocity]= initialize(L,N);
[dx,dy,r]=separation(N,L,x_position,y_position);
[Fx,Fy]=forceupdate(N,Fc,dx,dy,r);
%% Creating empty matrices for position, velocity, force and energies
x_instant=zeros(N,timesteps);
y_instant=zeros(N,timesteps);

vx_instant=zeros(N,timesteps);
vy_instant=zeros(N,timesteps);

Fx_instant=zeros(N,timesteps);
Fy_instant=zeros(N,timesteps);

Kinetic_Energy=zeros(1,timesteps);
Potential_Energy=zeros(1,timesteps);
Total_Energy=zeros(1,timesteps);

%% Saving the initial values of force and velocity
Fx_instant(:,1)=Fx;
Fy_instant(:,1)=Fy;

vx_instant(:,1)=x_velocity;
vy_instant(:,1)=y_velocity;


%% Main body

instant=1;%Setting up a counter for the timesteps

for i=1:length(t)
    %Save the current value of position as a new column
    x_instant(:,instant)=x_position;
    y_instant(:,instant)=y_position;
    for j=1:N
        %Update position and apply PBC
        x_position(j)=x_position(j)+vx_instant(j,i).*dt+(Fx_instant(j,i)/2)*dt2;
        if x_position(j)<0
        x_position(j)=x_position(j)+L;
        elseif x_position(j)>L
        x_position(j)=x_position(j)-L;
        end
        y_position(j)=y_position(j)+vy_instant(j,i).*dt+(Fy_instant(j,i)/2)*dt2;
        if y_position(j)<0
        y_position(j)=y_position(j)+L;
        elseif y_position(j)>L
        y_position(j)=y_position(j)-L;
        end
        %Calculate the new value of seperation for updated positions
        [dx,dy,r]=separation(N,L,x_position,y_position);
    end
    %Increment the timestep counter
    instant=instant+1;
    %Calcualte the new value of force and save it as a new column
    [Fx_instant(:,instant),Fy_instant(:,instant)]=forceupdate(N,Fc,dx,dy,r);
    %Calcualte the new values of velocities and save them as a new column
    [vx_instant(:,instant)]=velocityupdate(N,vx_instant(:,instant-1),dt,Fx_instant(:,(instant-1)),Fx_instant(:,instant));
    [vy_instant(:,instant)]=velocityupdate(N,vy_instant(:,instant-1),dt,Fy_instant(:,(instant-1)),Fy_instant(:,instant));
%     Ensuring that the average velocity remains 0
    vx_avg=sum(vx_instant(:,instant))/N;
    vy_avg=sum(vy_instant(:,instant))/N;
    for k = 1:N
    vx_instant(k,instant)=vx_instant(k,instant)-vx_avg;
    vy_instant(k,instant)=vy_instant(k,instant)-vy_avg;
    end
    v_magnitude=sqrt(vx_instant(:,instant).^2+vy_instant(:,instant).^2);
    %Calculate the value of energies at the current timestep
    [Kinetic_Energy(:,instant),Potential_Energy(:,instant),Total_Energy(:,instant)]=Energycalc(N,v_magnitude,r,rc,Fc,Vc);
end



%% Animation
f=figure("Position",[5,5,1080,540]);%Setting up the dimensions of the figure
frame=1;%creating a counter for the frame number
Colors=zeros(N,3);% Creating an empty matrix for the colormap
for i = 1:length(t)
    %We create a new frame only after 24 timesteps
    if mod(i,24)==0
    % extract the value of position and velocity at current timestep   
    x = x_instant(:,i);
    y = y_instant(:,i);
    v=sqrt(vx_instant(:,i).^2+vy_instant(:,i).^2);
    %Create a colormap according to the speeds of the particle
    for j= 1:N
    if v(j)<=max(v)/6
        Colors(j,:)=[0 0 1];
    elseif v(j)>max(v)/6 && v(j)<=2*max(v)/6
        Colors(j,:)=[0 1 0.8];
    elseif v(j)>2*max(v)/6 && v(j)<=3*max(v)/6
        Colors(j,:)=[0.4660 1 0.1880];
    elseif v(j)>3*max(v)/6 && v(j)<=4*max(v)/6
        Colors(j,:)=[0 1 0];
    elseif v(j)>4*max(v)/6 && v(j)<=5*max(v)/6
        Colors(j,:)=[1 0.5250 0.0980];
    elseif v(j)>5*max(v)/6 
        Colors(j,:)=[1 1 0];
    end
    end
    %Ploting the positions in the current frame as a scatter graph
    axis equal
    subplot(2,2,[1,3])
    s=scatter(x,y,65,Colors,'o','filled');
    colorbar('Ticks',[0,1],'TickLabels',{'Vmin','Vmax'});
    axis([0 L 0 L])
    set(gca,'color',[0 0 0])
    title('Main simulation box')
    xlabel('X Distance $$ \left(\frac{x}{\sigma}\right)$$','Interpreter','latex')
    ylabel('Y Distance $$ \left(\frac{y}{\sigma}\right)$$','Interpreter','latex')
    %plotting the normalised velocity distribution in the current frame as
    %a histogram
    subplot(2,2,2)
    histogram(v,6,"Normalization","pdf")
    title('Velocity probability distribution')
    xlabel('Velocity $$\left(v\cdot\left(\frac{m}{\epsilon}\right)^\frac{1}{2} \right)$$','Interpreter','latex')
    ylabel('Probability Density')
    %Plotting the energy values till the current frame
    subplot(2,2,4)
    plot(t(1:i),Total_Energy(1:i),'LineWidth',2,'Color','r');hold on;
    plot(t(1:i),Kinetic_Energy(1:i),'LineWidth',2,'Color','b');hold on
    plot(t(1:i),Potential_Energy(1:i),'LineWidth',2,'Color','g');
    axis ([0 tmax 0 0.5]);
    axis square
    legend('Total Energy','Kinetic Energy','Potential Energy','Location','bestoutside')
    title('Energy Evolution')
    xlabel('Time $$\left(t\cdot\left(\frac{\epsilon}{m\sigma^{2}}\right)^\frac{1}{2} \right)$$','Interpreter','latex')
    ylabel('Energy per particle $$\left(\frac{E}{N\cdot\epsilon}\right)$$','Interpreter','latex')
    %save the current frame in a flipbook
    flipbook(frame) = getframe(gcf);
    %Increment the frame value by 1
    frame=frame+1;
    end
end
%Save the frames as an mp4 animation
Fellini=VideoWriter('MD_N100_T10', 'MPEG-4');
Fellini.FrameRate=25;
Fellini.Quality=60;
open(Fellini);
writeVideo(Fellini,flipbook);
close(Fellini);
