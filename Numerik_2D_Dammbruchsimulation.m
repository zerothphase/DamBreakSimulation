clear;
clc;
close all;

%Inputs
N = 100;     
Nx = N;     %number of cells in x direction
Ny = N;     %number of cells in y direction
Nx_edges = Nx +1; %Number of cell edges
Ny_edges = Ny +1;
Lx1 = -1;   %range of domain in x direction
Lx2 = 1;
Ly1 = -1;   %range of domain in y direction
Ly2 = 1;
CFL = 0.8;
g = 10;
tmax = 0.5; %end time of the simulation

%Cell edges 
x_edges = linspace(Lx1,Lx2,Nx_edges);
y_edges = linspace(Ly1,Ly2,Ny_edges);
dx = x_edges(2)-x_edges(1);
dy = y_edges(2)-y_edges(1);

%Cell centers,Nx, Ny
x_centers = Lx1+dx/2 : dx : Lx2;
y_centers = Ly1+dy/2 : dy : Ly2;

[X,Y] = meshgrid(x_centers,y_centers);

domain = ones(Nx,Ny);   %Additional matrix to determine the computational domain.

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
h = zeros(Ny,Nx);
u = zeros(Ny,Nx);
v = zeros(Ny,Nx);


%Filter out regions between -0.2 < x < 0.2 and y > 0.2 or y < -0.2, set the
%element of "domain" to 0.
for i = 1:Nx
    for j = 1:Ny
        if X(j,i)>-0.2 && X(j,i)<0.2 && (Y(j,i)<-0.2 || Y(j,i)>0.2)
            domain(j,i) = 0;
        else
        end
    end
end

%Initialize
t = 0;
for i=1:Nx
    for j=1:Ny
        if X(j,i) < -0.2
            h(j,i) = 0.5;
            u(j,i) = 0;
            v(j,i) = 0;
        else
            h(j,i) = 1;
            u(j,i) = 0;
            v(j,i) = 0;
        end
    end
end

U1_new = h;
U2_new = h.*u;
U3_new = h.*v;
F1 = h.*u;
F2 = h.*u.^2+1/2*g*h.^2;
F3 = h.*u.*v;
G1 = h.*v;
G2 = h.*u.*v;
G3 = h.*v.^2+1/2*g*h.^2;

fig = figure;
filename = 'testAnimated.gif';
ind = 1;
while t < tmax
    %Compute global timestep
    dt = CFL/2*min(dx/max(max(abs(u)+sqrt(g*h))) , dy/max(max(abs(v)+sqrt(g*h))));
    if t+dt>tmax
        dt=tmax-t;
    end

    %horizontal fluxes
    for i = 1:Nx+1
        for j = 1:Ny
            if i == 1 %Outflow BC
                f1_edges(j,i) = 1/2*(F1(j,i)+F1(j,i))-dx/2/2/dt*CFL*(U1_new(j,i)-U1_new(j,i)); 
                f2_edges(j,i) = 1/2*(F2(j,i)+F2(j,i))-dx/2/2/dt*CFL*(U2_new(j,i)-U2_new(j,i));
                f3_edges(j,i) = 1/2*(F3(j,i)+F3(j,i))-dx/2/2/dt*CFL*(U3_new(j,i)-U3_new(j,i));
            elseif i == Nx+1 %Right Wall u_g = -u_r ->F1,F3,U2 -ve
                f1_edges(j,i) = 1/2*(F1(j,i-1)-F1(j,i-1))-dx/2/2/dt*CFL*(U1_new(j,i-1)-U1_new(j,i-1));
                f2_edges(j,i) = 1/2*(F2(j,i-1)+F2(j,i-1))-dx/2/2/dt*CFL*(-U2_new(j,i-1)-U2_new(j,i-1));
                f3_edges(j,i) = 1/2*(F3(j,i-1)-F3(j,i-1))-dx/2/2/dt*CFL*(U3_new(j,i-1)-U3_new(j,i-1));
            elseif domain(j,i-1) == 1 && domain(j,i) == 0 %Right Wall (x = -0.2)
                f1_edges(j,i) = 1/2*(F1(j,i-1)-F1(j,i-1))-dx/2/2/dt*CFL*(U1_new(j,i-1)-U1_new(j,i-1));
                f2_edges(j,i) = 1/2*(F2(j,i-1)+F2(j,i-1))-dx/2/2/dt*CFL*(-U2_new(j,i-1)-U2_new(j,i-1));
                f3_edges(j,i) = 1/2*(F3(j,i-1)-F3(j,i-1))-dx/2/2/dt*CFL*(U3_new(j,i-1)-U3_new(j,i-1));
            elseif domain(j,i-1) == 0 && domain(j,i) == 0
                f1_edges(j,i) = 0; 
                f2_edges(j,i) = 0;
                f3_edges(j,i) = 0;
            elseif domain(j,i-1) == 0 && domain(j,i) == 1 %Left Wall (x = 0.2)
                f1_edges(j,i) = 1/2*(-F1(j,i)+F1(j,i))-dx/2/2/dt*CFL*(U1_new(j,i)-U1_new(j,i));
                f2_edges(j,i) = 1/2*(F2(j,i)+F2(j,i))-dx/2/2/dt*CFL*(U2_new(j,i)+U2_new(j,i));
                f3_edges(j,i) = 1/2*(-F3(j,i)+F3(j,i))-dx/2/2/dt*CFL*(U3_new(j,i)-U3_new(j,i)); 
            else
                f1_edges(j,i) = 1/2*(F1(j,i-1)+F1(j,i))-dx/2/2/dt*CFL*(U1_new(j,i)-U1_new(j,i-1));
                f2_edges(j,i) = 1/2*(F2(j,i-1)+F2(j,i))-dx/2/2/dt*CFL*(U2_new(j,i)-U2_new(j,i-1));
                f3_edges(j,i) = 1/2*(F3(j,i-1)+F3(j,i))-dx/2/2/dt*CFL*(U3_new(j,i)-U3_new(j,i-1));           
            end
        end
    end

    for j = 1:Ny+1
        for i = 1:Nx
            if j == 1 && domain(j,i) == 1 %Wall
                g1_edges(j,i) = 1/2*(-G1(j,i)+G1(j,i))-dy/2/2/dt*CFL*(U1_new(j,i)-U1_new(j,i));
                g2_edges(j,i) = 1/2*(-G2(j,i)+G2(j,i))-dy/2/2/dt*CFL*(U2_new(j,i)-U2_new(j,i));
                g3_edges(j,i) = 1/2*(G3(j,i)+G3(j,i))-dy/2/2/dt*CFL*(U3_new(j,i)+U3_new(j,i));
            elseif j == 1 && domain(j,i) == 0 
                g1_edges(j,i) = 0;
                g2_edges(j,i) = 0;
                g3_edges(j,i) = 0;
            elseif j == Ny+1 && domain(j-1,i) == 1 %Wall
                g1_edges(j,i) = 1/2*(G1(j-1,i)-G1(j-1,i))-dy/2/2/dt*CFL*(U1_new(j-1,i)-U1_new(j-1,i));
                g2_edges(j,i) = 1/2*(G2(j-1,i)-G2(j-1,i))-dy/2/2/dt*CFL*(U2_new(j-1,i)-U2_new(j-1,i));
                g3_edges(j,i) = 1/2*(G3(j-1,i)+G3(j-1,i))-dy/2/2/dt*CFL*(-U3_new(j-1,i)-U3_new(j-1,i));
            elseif j == Ny+1 && domain(j-1,i) == 0
                g1_edges(j,i) = 0;
                g2_edges(j,i) = 0;
                g3_edges(j,i) = 0;
            elseif domain(j-1,i) == 1 && domain(j,i) == 0 %Wall (y = -0.2)
                g1_edges(j,i) = 1/2*(G1(j-1,i)-G1(j-1,i))-dy/2/2/dt*CFL*(U1_new(j-1,i)-U1_new(j-1,i));
                g2_edges(j,i) = 1/2*(G2(j-1,i)-G2(j-1,i))-dy/2/2/dt*CFL*(U2_new(j-1,i)-U2_new(j-1,i));
                g3_edges(j,i) = 1/2*(G3(j-1,i)+G3(j-1,i))-dy/2/2/dt*CFL*(-U3_new(j-1,i)-U3_new(j-1,i));
            elseif domain(j-1,i) == 0 && domain(j,i) == 0
                g1_edges(j,i) = 0;
                g2_edges(j,i) = 0;
                g3_edges(j,i) = 0;
            elseif domain(j-1,i) == 0 && domain(j,i) == 1 %Wall (y = 0.2)
                g1_edges(j,i) = 1/2*(-G1(j,i)+G1(j,i))-dy/2/2/dt*CFL*(U1_new(j,i)-U1_new(j,i));
                g2_edges(j,i) = 1/2*(-G2(j,i)+G2(j,i))-dy/2/2/dt*CFL*(U2_new(j,i)-U2_new(j,i));
                g3_edges(j,i) = 1/2*(G3(j,i)+G3(j,i))-dy/2/2/dt*CFL*(U3_new(j,i)+U3_new(j,i));
            else
                g1_edges(j,i) = 1/2*(G1(j-1,i)+G1(j,i))-dy/2/2/dt*CFL*(U1_new(j,i)-U1_new(j-1,i));
                g2_edges(j,i) = 1/2*(G2(j-1,i)+G2(j,i))-dy/2/2/dt*CFL*(U2_new(j,i)-U2_new(j-1,i));
                g3_edges(j,i) = 1/2*(G3(j-1,i)+G3(j,i))-dy/2/2/dt*CFL*(U3_new(j,i)-U3_new(j-1,i));
            end
        end
    end

        for i = 1:Nx
            for j = 1:Ny
                if domain(j,i) == 1
                    df1_edges(j,i) = f1_edges(j,i+1)-f1_edges(j,i);
                    df2_edges(j,i) = f2_edges(j,i+1)-f2_edges(j,i);
                    df3_edges(j,i) = f3_edges(j,i+1)-f3_edges(j,i);
                    dg1_edges(j,i) = g1_edges(j+1,i)-g1_edges(j,i);
                    dg2_edges(j,i) = g2_edges(j+1,i)-g2_edges(j,i);
                    dg3_edges(j,i) = g3_edges(j+1,i)-g3_edges(j,i);
                else
                    df1_edges(j,i) = 0;
                    df2_edges(j,i) = 0;
                    df3_edges(j,i) = 0;
                    dg1_edges(j,i) = 0;
                    dg2_edges(j,i) = 0;
                    dg3_edges(j,i) = 0;                
                end
            end
        end
    % Compute new U,h,u,v,F,G
    U1_new = U1_new - dt/dx*df1_edges - dt/dy*dg1_edges;
    U2_new = U2_new - dt/dx*df2_edges - dt/dy*dg2_edges;
    U3_new = U3_new - dt/dx*df3_edges - dt/dy*dg3_edges;

    h = U1_new;
    u = U2_new./U1_new;
    v = U3_new./U1_new;
    
    F1 = h.*u;
    F2 = h.*u.^2+1/2*g*h.^2;
    F3 = h.*u.*v;    
    G1 = h.*v;
    G2 = h.*u.*v;
    G3 = h.*v.^2+1/2*g*h.^2;
    t = t+dt;
    
    s = surf(X,Y,h,'EdgeColor','none');
    axis([-1 1 -1 1 0.4 1.2]);
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    zlabel('h','FontSize',14);
    text(1.1,1.1,1.1,['t = ',num2str(t,4),'s']);
    drawnow 
      % Capture the plot as an image 
      frame = getframe(fig); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if ind == 1
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',dt); 
      end 
    ind = ind + 1;
    %pause(dt);
end

for i=1:Nx
    for j=1:Ny
        if domain(j,i) == 0
            h(j,i) = NaN;
            u(j,i) = NaN;
            v(j,i) = NaN;
        else
        end
    end
end

%display maximum and minimum results
hmax = max(h(:));
umax = max(u(:));
vmax = max(v(:));
hmin = min(h(:));
umin = min(u(:));
vmin = min(v(:));

