function [rec_grid_conv_clip,theta_range2,r_grid,N] = cal_pattern(fc,N_subc)
c = physconst('Lightspeed');

% for fc_idx = 1:numel(fc_list)
%     fc = fc_list(fc_idx);
% fc = 10e9;


lambda = c/fc;
% N = 256;
N_UE = 1;
d = lambda/2;
% aperture = (N-1)*d; % ~1m
k = 2*pi/lambda;
% N_subc = 200;
%     for N_sub = [50,100,200]
f_SCS = 960e3;
BW = f_SCS*N_subc;

cylindrical_radius = 1;
R = cylindrical_radius;
alpha = pi/3;
theta_min = pi/2-alpha;
ref_source.x = 0;
ref_source.y = -R*cos(alpha); % these are virtual locations.
% theta_min = asin(-ref_source.y/cylindrical_radius);
% alpha = pi/2-theta_min;
N_min = ceil(4*alpha*R/lambda);
N = N_min;

R_min = R+1;
R_max = R+20;


R_UE = 0;
while R_UE<R_min
    R_UE = sqrt(rand())*R_max;
end
theta_UE = rand()*2*alpha*0.9+pi/2-alpha;
R_UE = 10;
theta_UE = pi/3;
UE.x = R_UE*cos(theta_UE);
UE.y = R_UE*sin(theta_UE)+ref_source.y;


%% coord setup
theta_range = linspace(pi/2-alpha, pi/2+alpha, N);
UCA.theta = theta_range';
UCA.x = R.*cos(UCA.theta);
UCA.y = R.*sin(UCA.theta)+ref_source.y;

recv_pattern = zeros(N,N_subc);
for subc=1:N_subc
    f_subc = fc+subc*f_SCS;
    lambda_subc = c/f_subc;
    k_subc = 2*pi/lambda_subc;
    recv_pattern(:,subc) = cal_recv_pattern(UCA.x,UCA.y,UE.x,UE.y,k_subc,1);
end
%% rec
G_d = 200;
delta_theta = 2*alpha/N;
theta_range2 = linspace(-(N-1),N-1,2*N-1)*delta_theta+delta_theta;
rec_grid_conv_k = zeros(3*N-2,G_d,N_subc);


r_grid = linspace(R_min,R_max,G_d);%cylindrical_radius;
for subc=1:N_subc
    f_subc = fc+subc*f_SCS;
    lambda_subc = c/f_subc;
    k_subc = 2*pi/lambda_subc;
    for r_idx = 1:numel(r_grid)
        r_trg = r_grid(r_idx);
        conv_kernel_k = exp(1j*k_subc*sqrt(R^2+r_trg^2-2*R*r_trg*cos(theta_range2)));
        res_k = conv(recv_pattern(:,subc),conv_kernel_k);
        rec_grid_conv_k(:,r_idx,subc) = res_k;
    end
end
rec_grid_conv_clip = rec_grid_conv_k(N:end-N+1,:,:);

end

