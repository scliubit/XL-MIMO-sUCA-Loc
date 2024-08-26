clear all
close all
%% ULA and UCA detection
%% ULA part
scatterer.phi = [pi/3,pi/12*7];
scatterer.r = [4,7];
scatterer.x = scatterer.r.*cos(scatterer.phi);
scatterer.y = scatterer.r.*sin(scatterer.phi);

fc = 3.5e9;
N = 64;
c = physconst('Lightspeed');
lambda = c/fc;
aperture = (N-1)*lambda/2;
kappa = 2*pi/lambda;

N_grid = 400;
x_grid = linspace(-4,4,N_grid);
y_grid = linspace(2,10,N_grid);

ULA.x = linspace(-aperture/2,aperture/2,N)';
ULA.y = zeros(N,1);

recv_ULA = sum(cal_recv_pattern(ULA.x,ULA.y,scatterer.x,scatterer.y,kappa,1),2);
rec_grid = zeros(numel(x_grid),numel(y_grid));
for x_idx = 1:numel(x_grid)
    for y_idx = 1:numel(y_grid)
        x_i = x_grid(x_idx);
        y_i = y_grid(y_idx);
        distance = sqrt((x_i-ULA.x).^2+(y_i-ULA.y).^2);
        rec_grid(x_idx,y_idx)=sum(recv_ULA.*exp(1j*kappa*distance));
    end
end
%% UCA part
alpha = pi/3;
theta_min = pi/2-alpha;
theta_max = pi/2+alpha;

c = physconst('Lightspeed');
lambda = c/fc;
N_UE = 1;
d = lambda/2;
k = 2*pi/lambda;
cylindrical_radius = 1;
R = cylindrical_radius;

UCA.theta = linspace(theta_min,theta_max,N)';
UCA.r = R+zeros(N,1);
UCA.x = (UCA.r).*cos(UCA.theta);
UCA.y = (UCA.r).*sin(UCA.theta)-0.5;

r_grid = linspace(2,10,N_grid);
theta_grid = linspace(theta_min,theta_max,N_grid);
recv_UCA = sum(cal_recv_pattern(UCA.x,UCA.y,scatterer.x,scatterer.y,kappa,1),2);

rec_grid2 = zeros(numel(r_grid),numel(theta_grid));
for r_idx = 1:numel(r_grid)
    for t_idx = 1:numel(theta_grid)
        r_i = r_grid(r_idx);
        t_i = theta_grid(t_idx);
        distance = sqrt(R^2+r_i^2-2*R*r_i*cos(t_i-UCA.theta));
        rec_grid2(r_idx,t_idx)=sum(recv_UCA.*exp(1j*kappa*distance));
    end
end

figure
subplot(221)
surf(x_grid,y_grid,abs(rec_grid)')
shading interp

xlabel('$x$ (m)',Interpreter='latex')
ylabel('$y$ (m)',Interpreter='latex')
title('(a) ULA Ambiguity Function',Interpreter='latex')

set(gca,'fontsize',14);
set(gca,'fontname','Times New Roman')

view([0 0 1])
axis('tight')
box on

subplot(223)

rec_grid_sum = sum(abs(rec_grid),2);

plot(x_grid, flipud(rec_grid_sum)/max(rec_grid_sum(:)))
xlabel('$x$ (m)',Interpreter='latex')
ylabel('Normalized Amplitude',Interpreter='latex')
title('(c) $x$-axis Projection',Interpreter='latex')

set(gca,'fontsize',14);
set(gca,'fontname','Times New Roman')
grid on
box on

subplot(222)
surf(theta_grid,r_grid,abs(rec_grid2))
shading interp

xlabel('$\varphi$ (rad)',Interpreter='latex')
ylabel('$r$ (m)',Interpreter='latex')
title('(b) sUCA Ambiguity Function',Interpreter='latex')

set(gca,'fontsize',14);
set(gca,'fontname','Times New Roman')

view([0 0 1])
axis('tight')
box on

subplot(224)
rec_grid2_sum = sum(abs(rec_grid2),1);

plot(theta_grid, flipud(rec_grid2_sum)/max(rec_grid2_sum(:)))
xlabel('$\varphi$ (rad)',Interpreter='latex')
ylabel('Normalized Amplitude',Interpreter='latex')
title('(d) $\varphi$-axis Projection',Interpreter='latex')

set(gca,'fontsize',14);
set(gca,'fontname','Times New Roman')
grid on
box on
axis('tight')

% set(gcf,'renderer','Painters')
% plot(rec_grid2_sum)


