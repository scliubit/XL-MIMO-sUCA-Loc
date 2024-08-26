clear all
close all
[rec_grid_conv_clip,theta_range2,r_grid,N] = cal_pattern(10e9,100);
%% const table
R = 1;
R_min = R+1;
R_max = R+20;
alpha = pi/3;
G_d = 200;
ref_source.x = 0;
ref_source.y = -R*cos(alpha); % these are virtual locations.
theta_range = linspace(pi/2-alpha, pi/2+alpha, N);
UCA.theta = theta_range';
UCA.x = R.*cos(UCA.theta);
UCA.y = R.*sin(UCA.theta)+ref_source.y;
R_UE = 10;
theta_UE = pi/3;
UE.x = R_UE*cos(theta_UE);
UE.y = R_UE*sin(theta_UE)+ref_source.y;


angle_grid2 = linspace(pi/2-alpha,pi/2+alpha,N);
rec_grid_conv_clip_sum = sum(rec_grid_conv_clip,3);
rec_grid_conv_clip_norm = rec_grid_conv_clip_sum/max(abs(rec_grid_conv_clip_sum(:)));
%% detect

% rec_grid_conv_clip_sum = abs(sum(rec_grid_conv_clip,3));
[~,idx]=max(abs(rec_grid_conv_clip_norm(:)));
[idx1,idx2] = ind2sub(size(rec_grid_conv_clip_norm),idx);
angle_est = angle_grid2(idx1);
r_est = r_grid(idx2);%-ref_source.y;
% disp()
%% plot polar coord on cartesian coord

figure
subplot(121)
hold on
angle_grid = linspace(pi/2-alpha,pi/2+alpha,N);
r_grid = linspace(R_min,R_max,G_d);
[X,Y] = meshgrid(r_grid,angle_grid);
for ang_idx = 1:numel(angle_grid2)
    for r_idx = 1:numel(r_grid)
        X(ang_idx,r_idx) = r_grid(r_idx)*cos(angle_grid2(ang_idx)-pi/2)+ref_source.y;
        Y(ang_idx,r_idx) = r_grid(r_idx)*sin(angle_grid2(ang_idx)-pi/2);
    end
end
Z = abs(rec_grid_conv_clip_norm);
surf(X,Y,Z);

shading interp
xlabel('y (m)',Interpreter='latex')
ylabel('x (m)',Interpreter='latex')
zlabel('Normalized Amplitude',Interpreter='latex')
title('(a) Reconstruction with sUCA',Interpreter='latex')
view(0,90)
colorbar()
% axis('tight')

grid off
box on
plot3(R*cos(UCA.theta-pi/2)+ref_source.y,R*sin(UCA.theta-pi/2),ones(numel(UCA.theta),1),'g-',LineWidth=2)
plot3(R_UE*cos(theta_UE-pi/2)+ref_source.y,R_UE*sin(theta_UE-pi/2),1,'ro',MarkerSize=15,LineWidth=2)
legend({'','sUCA','UE $(r_0, \varphi_0)$'},Interpreter="latex",Location="northeast")
set(gca,'fontsize',14)
set(gca,'fontname','Times New Roman')
axis('equal')

%% upper right

subplot(222)
hold on
% 3.5G
[rec_grid_conv_clip,~,r_grid,N] = cal_pattern(3.5e9,100);
angle_grid = linspace(pi/2-alpha,pi/2+alpha,N);
rec_grid_conv_clip_sum = sum(rec_grid_conv_clip,3);
rec_grid_conv_clip_norm = rec_grid_conv_clip_sum/max(abs(rec_grid_conv_clip_sum(:)));
[~,idx]=max(abs(rec_grid_conv_clip_norm(:)));
[idx1,idx2] = ind2sub(size(rec_grid_conv_clip_norm),idx);
angle_est = angle_grid2(idx1);
r_est = r_grid(idx2);%-ref_source.y;
plot(angle_grid,abs(rec_grid_conv_clip_norm(:,idx2)),'-.',LineWidth=2)
% 10G
[rec_grid_conv_clip,~,r_grid,N] = cal_pattern(10e9,100);
angle_grid = linspace(pi/2-alpha,pi/2+alpha,N);
rec_grid_conv_clip_sum = sum(rec_grid_conv_clip,3);
rec_grid_conv_clip_norm = rec_grid_conv_clip_sum/max(abs(rec_grid_conv_clip_sum(:)));
[~,idx]=max(abs(rec_grid_conv_clip_norm(:)));
[idx1,idx2] = ind2sub(size(rec_grid_conv_clip_norm),idx);
angle_est = angle_grid2(idx1);
r_est = r_grid(idx2);%-ref_source.y;
plot(angle_grid,abs(rec_grid_conv_clip_norm(:,idx2)),'-.',LineWidth=2,Color="#EDB120")
% 28G
[rec_grid_conv_clip,~,r_grid,N] = cal_pattern(28e9,100);
angle_grid = linspace(pi/2-alpha,pi/2+alpha,N);
rec_grid_conv_clip_sum = sum(rec_grid_conv_clip,3);
rec_grid_conv_clip_norm = rec_grid_conv_clip_sum/max(abs(rec_grid_conv_clip_sum(:)));
[~,idx]=max(abs(rec_grid_conv_clip_norm(:)));
[idx1,idx2] = ind2sub(size(rec_grid_conv_clip_norm),idx);
angle_est = angle_grid2(idx1);
r_est = r_grid(idx2);%-ref_source.y;
plot(angle_grid,abs(rec_grid_conv_clip_norm(:,idx2)),'-.',LineWidth=2,Color="#7E2F8E")
legend({'$f_c=3.5\,$GHz','$f_c=10\,$GHz','$f_c=28\,$GHz'},Interpreter="latex",Location="northeast")
axis('tight')
xlabel('$\varphi$ (rad)',Interpreter='latex')
ylabel('Normalized Amplitude',Interpreter='latex')
title('(b) Angular Spectrum $r = r_0$',Interpreter='latex')
set(gca,'fontsize',14)
set(gca,'fontname','Times New Roman')
grid on
box on


subplot(224)
hold on
% angular
% BW 24M
[rec_grid_conv_clip,~,r_grid,N] = cal_pattern(10e9,100);
angle_grid = linspace(pi/2-alpha,pi/2+alpha,N);
rec_grid_conv_clip_sum = sum(rec_grid_conv_clip,3);
rec_grid_conv_clip_norm = rec_grid_conv_clip_sum/max(abs(rec_grid_conv_clip_sum(:)));
[~,idx]=max(abs(rec_grid_conv_clip_norm(:)));
[idx1,idx2] = ind2sub(size(rec_grid_conv_clip_norm),idx);
angle_est = angle_grid2(idx1);
r_est = r_grid(idx2);%-ref_source.y;
plot(r_grid,abs(rec_grid_conv_clip_norm(idx1,:)),'-',LineWidth=2,Color="#FF0000")
% BW 48M
[rec_grid_conv_clip,~,r_grid,N] = cal_pattern(10e9,150);
angle_grid = linspace(pi/2-alpha,pi/2+alpha,N);
rec_grid_conv_clip_sum = sum(rec_grid_conv_clip,3);
rec_grid_conv_clip_norm = rec_grid_conv_clip_sum/max(abs(rec_grid_conv_clip_sum(:)));
[~,idx]=max(abs(rec_grid_conv_clip_norm(:)));
[idx1,idx2] = ind2sub(size(rec_grid_conv_clip_norm),idx);
angle_est = angle_grid2(idx1);
r_est = r_grid(idx2);%-ref_source.y;
plot(r_grid,abs(rec_grid_conv_clip_norm(idx1,:)),'-',LineWidth=2,Color="#0072BD")
% BW 96M
[rec_grid_conv_clip,~,r_grid,N] = cal_pattern(10e9,200);
angle_grid = linspace(pi/2-alpha,pi/2+alpha,N);
rec_grid_conv_clip_sum = sum(rec_grid_conv_clip,3);
rec_grid_conv_clip_norm = rec_grid_conv_clip_sum/max(abs(rec_grid_conv_clip_sum(:)));
[~,idx]=max(abs(rec_grid_conv_clip_norm(:)));
[idx1,idx2] = ind2sub(size(rec_grid_conv_clip_norm),idx);
angle_est = angle_grid2(idx1);
r_est = r_grid(idx2);%-ref_source.y;
plot(r_grid,abs(rec_grid_conv_clip_norm(idx1,:)),'-',LineWidth=2,Color="#EDB120")


axis('tight')
xlabel('$r$ (m)',Interpreter='latex')
ylabel('Normalized Amplitude',Interpreter='latex')
legend({'$K=100$','$K=150$','$K=200$'},Interpreter="latex",Location="northeast")
title('(c) Distance Spectrum $\varphi = \varphi_0$',Interpreter='latex')
set(gca,'fontsize',14)
set(gca,'fontname','Times New Roman')
grid on
box on


% set(cf,'renderer','Painters')



