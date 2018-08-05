clear;
clc;
close all;

%Inputs
N = 25;     
Nx = N;     %number of cells in x direction
Ny = N;     %number of cells in y direction
Lx1 = -1;   %range of domain in x direction
Lx2 = 1;
Ly1 = -1;   %range of domain in y direction
Ly2 = 1;
CFL = 0.8;
g = 10;
tmax = 1;   %end time of the simulation

Nx_edges = Nx +1; %Number of cell edges
Ny_edges = Ny +1;

%Define manufactured solutions
h_ex = @(x,y,t) (cos(pi*(x+y+t))+2);
u_ex = @(x,y,t) cos(pi*(x+y+t));
v_ex = @(x,y,t) cos(pi*(x+y+t));

f1_edges = zeros(Ny,Nx+1);
f2_edges = zeros(Ny,Nx+1);
f3_edges = zeros(Ny,Nx+1);
g1_edges = zeros(Ny+1,Nx);
g2_edges = zeros(Ny+1,Nx);
g3_edges = zeros(Ny+1,Nx);
df1_edges = zeros(Ny,Nx);
df2_edges = zeros(Ny,Nx);
df3_edges = zeros(Ny,Nx);
dg1_edges = zeros(Ny,Nx);
dg2_edges = zeros(Ny,Nx);
dg3_edges = zeros(Ny,Nx);

%Cell edges 
x_edges = linspace(Lx1,Lx2,Nx_edges);
y_edges = linspace(Ly1,Ly2,Ny_edges);
dx = x_edges(2)-x_edges(1);
dy = y_edges(2)-y_edges(1);

%Cell centers,Nx, Ny
x_centers = Lx1+dx/2 : dx : Lx2;
y_centers = Ly1+dy/2 : dy : Ly2;

[X,Y] = meshgrid(x_centers,y_centers);

%Initialize
t = 0;
h = h_ex(X,Y,0);
u = u_ex(X,Y,0);
v = v_ex(X,Y,0);
surf(X,Y,h);
U1_new = h;
U2_new = h.*u;
U3_new = h.*v;
F1 = h.*u;
F2 = h.*u.^2+1/2*g*h.^2;
F3 = h.*u.*v;
G1 = h.*v;
G2 = h.*u.*v;
G3 = h.*v.^2+1/2*g*h.^2;
Sa = -pi*sin(pi*(X+Y+t)).*(4*cos(pi*(X+Y+t))+5);
Sb = -pi*sin(pi*(X+Y+t)).*((10+g)*cos(pi*(X+Y+t))+6*(cos(pi*(X+Y+t))).^2+2+2*g);
Sc = Sb;


%Transient simulation

while t < tmax
    %Compute global timestep dt_2D
    dt = CFL/2*min(dx/max(max(abs(u)+sqrt(g*h))) , dy/max(max(abs(v)+sqrt(g*h))));
    if t+dt>tmax
        dt=tmax-t;
    end
    
    %horizontal fluxes
    for i = 1:Nx_edges
        for j = 1:Ny
            %left periodic BC
            if i == 1
                f1_edges(j,i) = 1/2*(F1(j,Nx)+F1(j,i))-dx/2/2/dt*CFL*(U1_new(j,i)-U1_new(j,Nx)); %1D timestep should be used to compute the numerical flux: dt_1D = dt_2D
                f2_edges(j,i) = 1/2*(F2(j,Nx)+F2(j,i))-dx/2/2/dt*CFL*(U2_new(j,i)-U2_new(j,Nx));
                f3_edges(j,i) = 1/2*(F3(j,Nx)+F3(j,i))-dx/2/2/dt*CFL*(U3_new(j,i)-U3_new(j,Nx));
            %right periodic BC
            elseif i == Nx_edges
                f1_edges(j,i) = 1/2*(F1(j,i-1)+F1(j,1))-dx/2/2/dt*CFL*(U1_new(j,1)-U1_new(j,i-1)); 
                f2_edges(j,i) = 1/2*(F2(j,i-1)+F2(j,1))-dx/2/2/dt*CFL*(U2_new(j,1)-U2_new(j,i-1));
                f3_edges(j,i) = 1/2*(F3(j,i-1)+F3(j,1))-dx/2/2/dt*CFL*(U3_new(j,1)-U3_new(j,i-1));
            %Lax-Friedrich-Fluss
            else
                f1_edges(j,i) = 1/2*(F1(j,i-1)+F1(j,i))-dx/2/2/dt*CFL*(U1_new(j,i)-U1_new(j,i-1)); 
                f2_edges(j,i) = 1/2*(F2(j,i-1)+F2(j,i))-dx/2/2/dt*CFL*(U2_new(j,i)-U2_new(j,i-1));
                f3_edges(j,i) = 1/2*(F3(j,i-1)+F3(j,i))-dx/2/2/dt*CFL*(U3_new(j,i)-U3_new(j,i-1));
            end
        end
    end

    %vertical fluxes
    for j = 1:Ny_edges
        for i = 1:Nx
            %top periodic BC
            if j == 1
                g1_edges(j,i) = 1/2*(G1(Ny,i)+G1(j,i))-dy/2/2/dt*CFL*(U1_new(j,i)-U1_new(Ny,i)); %Ny+1*Nx
                g2_edges(j,i) = 1/2*(G2(Ny,i)+G2(j,i))-dy/2/2/dt*CFL*(U2_new(j,i)-U2_new(Ny,i));
                g3_edges(j,i) = 1/2*(G3(Ny,i)+G3(j,i))-dy/2/2/dt*CFL*(U3_new(j,i)-U3_new(Ny,i));
            %bottom periodic BC
            elseif j == Ny_edges
                g1_edges(j,i) = 1/2*(G1(j-1,i)+G1(1,i))-dy/2/2/dt*CFL*(U1_new(1,i)-U1_new(j-1,i)); %Ny+1*Nx
                g2_edges(j,i) = 1/2*(G2(j-1,i)+G2(1,i))-dy/2/2/dt*CFL*(U2_new(1,i)-U2_new(j-1,i));
                g3_edges(j,i) = 1/2*(G3(j-1,i)+G3(1,i))-dy/2/2/dt*CFL*(U3_new(1,i)-U3_new(j-1,i));
            %Lax-Friedrich-Fluss
            else
                g1_edges(j,i) = 1/2*(G1(j-1,i)+G1(j,i))-dy/2/2/dt*CFL*(U1_new(j,i)-U1_new(j-1,i)); %Ny+1*Nx
                g2_edges(j,i) = 1/2*(G2(j-1,i)+G2(j,i))-dy/2/2/dt*CFL*(U2_new(j,i)-U2_new(j-1,i));
                g3_edges(j,i) = 1/2*(G3(j-1,i)+G3(j,i))-dy/2/2/dt*CFL*(U3_new(j,i)-U3_new(j-1,i));
            end
        end
    end

    for i = 1:Nx
        for j = 1:Ny
            df1_edges(j,i) = f1_edges(j,i+1)-f1_edges(j,i);
            df2_edges(j,i) = f2_edges(j,i+1)-f2_edges(j,i);
            df3_edges(j,i) = f3_edges(j,i+1)-f3_edges(j,i);
            dg1_edges(j,i) = g1_edges(j+1,i)-g1_edges(j,i);
            dg2_edges(j,i) = g2_edges(j+1,i)-g2_edges(j,i);
            dg3_edges(j,i) = g3_edges(j+1,i)-g3_edges(j,i);
        end
    end


    % Compute new U,h,u,v,F,G
    U1_new = U1_new - dt/dx*df1_edges - dt/dy*dg1_edges+ Sa*dt;
    U2_new = U2_new - dt/dx*df2_edges - dt/dy*dg2_edges+ Sb*dt;
    U3_new = U3_new - dt/dx*df3_edges - dt/dy*dg3_edges+ Sc*dt;

    h = U1_new;
    u = U2_new./U1_new;
    v = U3_new./U1_new;
    
    F1 = h.*u;
    F2 = h.*u.^2+1/2*g*h.^2;
    F3 = h.*u.*v;    
    G1 = h.*v;
    G2 = h.*u.*v;
    G3 = h.*v.^2+1/2*g*h.^2;
    Sa = -pi*sin(pi*(X+Y+t)).*(4*cos(pi*(X+Y+t))+5);
    Sb = -pi*sin(pi*(X+Y+t)).*((10+g)*cos(pi*(X+Y+t))+6*(cos(pi*(X+Y+t))).^2+2+2*g);
    Sc = Sb;
    
    
    
    % Plot
    surf(X,Y,h);
    axis([-1 1 -1 1 0.5 3.5]);
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    zlabel('h','FontSize',14);
    
    t = t+dt;
    pause(dt);
end

%display max error and rms error
    quaderror = (h-h_ex(X,Y,t)).^2;
    RMS = sqrt(mean(quaderror(:)))
    error = abs(h-h_ex(X,Y,t));
    MAX = max(error(:))

