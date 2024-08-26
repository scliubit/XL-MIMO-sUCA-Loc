clear all
close all
%% calculate incomplete bessel function


%% practical setup
fc_list = [3.5e9, 10e9, 28e9];
alpha = pi/3;
theta_min = pi/2-alpha;
theta_max = pi/2+alpha;
phi_list = theta_min:0.001:theta_max;

f_phi = zeros(numel(phi_list),numel(fc_list));
f_phi_discrete = zeros(numel(phi_list),numel(fc_list));
for fc_idx = 1:numel(fc_list)
    fc = fc_list(fc_idx);

    c = physconst('Lightspeed');
    lambda = c/fc;
    N_UE = 1;
    d = lambda/2;
    k = 2*pi/lambda;
    cylindrical_radius = 1;
    R = cylindrical_radius;

    ref_source.x = 0;
    ref_source.y = -R*cos(alpha); % these are virtual locations.
    N_min = ceil(4*alpha*R/lambda);
    N = N_min;

    theta_range = linspace(pi/2-alpha, pi/2+alpha, N);
    UCA.theta = theta_range';
    UCA.x = R.*cos(UCA.theta);
    UCA.y = R.*sin(UCA.theta)+ref_source.y;
    R_UE = 10;
    theta_UE = pi/3;
    UE.x = R_UE*cos(theta_UE);
    UE.y = R_UE*sin(theta_UE)+ref_source.y;



    R_min = R+1;
    R_max = R+20;
    theta_min = pi/2-alpha;
    theta_max = pi/2+alpha;

    %%
    % alpha = pi/3;

    scatterer.r = 5;
    scatterer.phi = pi/3;

    integral_core_trial = @(phi,theta) exp(-2j*k*cylindrical_radius*sin(theta-(phi+scatterer.phi)/2)*sin((phi-scatterer.phi)/2));
    integral_core_trial2 = @(phi,theta) exp(-2j*k*cylindrical_radius*sin(theta-(phi+scatterer.phi)/2)*sin((phi-scatterer.phi)/2))./(scatterer.r-cylindrical_radius*cos(theta-scatterer.phi));
    cnt = 1;

    for phi_i = phi_list
        f_phi(cnt,fc_idx) = integral(@(theta) integral_core_trial2(phi_i,theta),theta_min,theta_max);
        f_phi_discrete(cnt,fc_idx) = sum(exp(-2j*k*cylindrical_radius*sin(UCA.theta-(phi_i+scatterer.phi)/2)*sin((phi_i-scatterer.phi)/2))./(scatterer.r-cylindrical_radius*cos(UCA.theta-scatterer.phi)));
        cnt=cnt+1;
    end
end
figure
subplot(211)
lambda1 = c/fc_list(1);
t = (-alpha:0.001:alpha)-scatterer.phi;
% rangee = max(sin(t))-min(sin(t));
dphi = lambda1/R/2/sin(alpha);
dphi_grids = dphi/0.001;
[~,max_idx] = max(abs(f_phi(:,1)));
max_angle = phi_list(max_idx);
hold on
plot(phi_list,(abs(f_phi(:,1))/max(abs(f_phi(:,1)))),LineWidth=1)
plot(phi_list,(abs(f_phi_discrete(:,1))/max(abs(f_phi_discrete(:,1)))),'--',LineWidth=1)
plot(phi_list,(abs(f_phi(:,3))/max(abs(f_phi(:,3)))),LineWidth=1)
plot(phi_list,(abs(f_phi_discrete(:,3))/max(abs(f_phi_discrete(:,3)))),'--',LineWidth=1)
plot([max_angle-dphi,max_angle-dphi,max_angle+dphi,max_angle+dphi],[0,1,1,0],'r-.',LineWidth=1.5)
lambda2 = c/fc_list(3);
dphi = lambda2/R/sin(alpha)/2;
dphi_grids = dphi/0.001;
[~,max_idx] = max(abs(f_phi(:,3)));
max_angle = phi_list(max_idx);
plot([max_angle-dphi,max_angle-dphi,max_angle+dphi,max_angle+dphi],[0,1,1,0],'b-.',LineWidth=1.5)

grid on
box on
xlabel('$\varphi$ (rad)',Interpreter='latex')
ylabel('Normalized Amplitude',Interpreter='latex')
title("\bf (a) Angular Pattern at $r=r_\ell$",Interpreter="latex")
legend({'$f_c=3.5\,$GHz, Integral','$f_c=3.5\,$GHz, Discrete Sum', '$f_c=28\,$GHz, Integral','$f_c=28\,$GHz, Discrete Sum', '$3.5\,$GHz Main Lobe (Est.)','$28\,$GHz Lobe Width (Est.)'},Interpreter="latex",Location="northeast")
set(gca,'fontsize',14);
set(gca,'fontname','Times New Roman')
% axis('tight')
axis([max_angle-pi/40 max_angle+pi/15 0 1])

%% ========================= distance ====================
fc_list = [3.5e9, 28e9];
alpha = pi/3;
theta_min = pi/2-alpha;
theta_max = pi/2+alpha;
phi_list = theta_min:0.001:theta_max;
f_SCS = 480e3;
G_d = 200;
k_idx = 0;
k_list = [50,100,200];
f_r = zeros(G_d,numel(fc_list),numel(k_list),200);
f_r_discrete = zeros(G_d,numel(fc_list),numel(k_list),200);
for K = k_list
    k_idx = k_idx+1;
    
    for fc_idx = 1:numel(fc_list)
        fc = fc_list(fc_idx);

        c = physconst('Lightspeed');
        lambda = c/fc;
        N_UE = 1;
        d = lambda/2;
        k = 2*pi/lambda;
        cylindrical_radius = 1;
        R = cylindrical_radius;

        ref_source.x = 0;
        ref_source.y = -R*cos(alpha); % these are virtual locations.
        N_min = ceil(4*alpha*R/lambda);
        N = N_min;

        theta_range = linspace(pi/2-alpha, pi/2+alpha, N);
        UCA.theta = theta_range';
        UCA.x = R.*cos(UCA.theta);
        UCA.y = R.*sin(UCA.theta)+ref_source.y;
        R_UE = 10;
        theta_UE = pi/3;
        UE.x = R_UE*cos(theta_UE);
        UE.y = R_UE*sin(theta_UE)+ref_source.y;

        R_min = R+1;
        R_max = R+20;


        theta_min = pi/2-alpha;
        theta_max = pi/2+alpha;

        %%
        % alpha = pi/3;

        scatterer.r = 10;
        scatterer.phi = pi/3;

%         integral_core_trial = @(phi,theta) exp(-2j*k*cylindrical_radius*sin(theta-(phi+scatterer.phi)/2)*sin((phi-scatterer.phi)/2));
%         integral_core_trial2 = @(phi,theta) exp(-2j*k*cylindrical_radius*sin(theta-(phi+scatterer.phi)/2)*sin((phi-scatterer.phi)/2))./(scatterer.r-cylindrical_radius*cos(theta-scatterer.phi));
        cnt = 1;
        r_gird = linspace(R_min,R_max,G_d)';
        for kk = 1:K
            f_r(:,fc_idx,k_idx,kk) = f_r(:,fc_idx,k_idx,kk) + exp(1j*2*pi/c*(fc+kk*f_SCS).*(scatterer.r+R^2/scatterer.r/2-(r_gird+R^2./r_gird)));
        end
    end
end

f_r11 = abs(sum(f_r(:,1,1,:),4));
f_r11 = f_r11/max(f_r11(:));
f_r12 = abs(sum(f_r(:,1,2,:),4));
f_r12 = f_r12/max(f_r12(:));
f_r13 = abs(sum(f_r(:,1,3,:),4));
f_r13 = f_r13/max(f_r13(:));


subplot(212)
hold on
plot(r_gird,f_r11,'-',LineWidth=2)
plot(r_gird,f_r12,'-.',LineWidth=2)
plot(r_gird,f_r13,'--',LineWidth=2)
grid on
xlabel('$r$ (m)',Interpreter='latex')
ylabel('Normalized Amplitude',Interpreter='latex')
title("\bf (b) Distance Pattern at $\varphi=\varphi_\ell$",Interpreter="latex")

width2 = 2*c/(k_list(3)*f_SCS);
plot([scatterer.r-width2/2 scatterer.r-width2/2],[0 1],'--k','Linewidth',1.4)
plot([scatterer.r+width2/2 scatterer.r+width2/2],[0 1],'--k','Linewidth',1.4)
legend({'$K=50$','$K=100$','$K=200$'},Interpreter="latex",Location="northeast")



% plot(r_gird,f_r22)
set(gca,'fontsize',14);
set(gca,'fontname','Times New Roman')
axis([R_min R_max 0 1])
x = [(scatterer.r+width2/2)/(R_max-R_min)-0.07 (scatterer.r+width2/2)/(R_max-R_min)]-0.11;
y = [0.5 0.5]-0.2;
Xadj = 0.91;
annotation('textarrow',x,y,'String',{'$\rm Width$   ','$K=200$'},'FontSize',16,'Linewidth',2,FontName="times new roman",Interpreter="latex")
annotation('textarrow',-x+Xadj,y,'String','','FontSize',14,'Linewidth',2)
box on
