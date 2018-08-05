clear;
clc;
close all;

%Change case here to plot different results
%1: plot h, 2: plot hu
Case = 1;


%plot exact solution
load u_ex.dat;
x_exact = u_ex(:,1);
h_exact = u_ex(:,2);
hu_exact = u_ex(:,3);

if Case == 1
    plot(x_exact,h_exact,'-');
    xlabel('x','FontSize',14);
    ylabel('h','FontSize',14);
else
    plot(x_exact,hu_exact,'-');
    xlabel('x','FontSize',14);
    ylabel('hu','FontSize',14);
end

hold all;
%Grid
N = 25; %number of cells
%Start first simulation with N = 25, then repeat with N = 2*N for 2 times
for l = 1:3
    N_edges = N+1; %Number of cell edges
    Lx1 = -1;
    Lx2 = 1;

    %Cell edges (x_i-1/2),(x_i+1/2),...
    x_edges = linspace(Lx1,Lx2,N_edges);
    dx = x_edges(2)-x_edges(1);

    %Cell centers (x_i)(x_i+1),...
    x_centers = Lx1+dx/2 : dx : Lx2;

    %Initialize
    U1 = zeros(1,N);
    U2 = zeros(1,N);
    f_edges = zeros(1,N_edges);
    g_edges = zeros(1,N_edges);

    for i=1:N
        if x_centers(i) <= 0
            U1(i) = 1;
        else
            U1(i) = 0.5;
        end
    end

    h = U1;
    u = U2./U1;

    %Transient simulation

    %Inputs
    CFL = 0.8;
    g = 10; 
    tmax = 0.2;
    t = 0;


    while t <= tmax
        % numerical flux at the edges
        dt = CFL*dx/max(abs(u)+sqrt(abs(g*h)));
        
        for i = 1:N_edges
            if i == 1
                f_edges(i) = h(i)*u(i);
                g_edges(i) = h(i)*u(i)^2 + 1/2*g*h(i)^2;

            elseif i == N_edges
                f_edges(i) = h(i-1)*u(i-1);
                g_edges(i) = h(i-1)*u(i-1)^2 + 1/2*g*h(i-1)^2;
            else
                u_left = h(i-1);
                u_right = h(i);
                v_left = h(i-1)*u(i-1);
                v_right = h(i)*u(i);
                f_left = h(i-1)*u(i-1);
                f_right = h(i)*u(i);
                g_left = h(i-1)*u(i-1)^2 + 1/2*g*h(i-1)^2;
                g_right = h(i)*u(i)^2 + 1/2*g*h(i)^2;
                f_edges(i) = 1/2*(f_left+f_right)-dx/2/dt*CFL*(u_right-u_left);
                g_edges(i) = 1/2*(g_left+g_right)-dx/2/dt*CFL*(v_right-v_left);
            end
        end
        
        for i = 1:N
            df(i) = f_edges(i+1)-f_edges(i);
            dg(i) = g_edges(i+1)-g_edges(i);
        end

        % Compute new U
        U1 = U1 -dt/dx*df;
        U2 = U2 -dt/dx*dg;

        % Extract new u and h
        h = U1;
        u = U2./U1;
        t = t+dt;

    %     pause(0.01);    
    end
    if Case == 1
        plot(x_centers,h,'-');
    else
        plot(x_centers,h.*u,'-');
    end
    N = N*2; %25,50,100 %Number of cells for the next loop
end
legend('exact','N=25','N=50','N=100')

