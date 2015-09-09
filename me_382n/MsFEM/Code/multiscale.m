% Multiscale FEM Project
%
% Miles Gray
% Dylan Pederson
clear all; 

%% relevant physical properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% material that the medium is composed of
conductivity_material = 2.4;    % W/m-K
density_material = 2.2/1000*100^3;         % kg/m^3
specific_heat_material = 1.0*1000;   % J/kg-K;
% filler material (porosity component) 
conductivity_filler = 0.02;     % W/m-K
density_filler = 1.0;           % kg/m^3
specific_heat_filler = 1.0*1000;     % J/kg-K;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define a coarse grid
%%%%%%%%%%%%%%%%%%%%%
x_min = 0; % m
x_max = 1e-2; % m
spacing_coarse = 1e-4; % m
%%%%%%%%%%%%%%%%%%%%%
x_coarse = x_min:spacing_coarse:x_max;
coarse_npts = length(x_coarse);

% define a fine grid
%%%%%%%%%%%%%%%%%%%%%
spacing_fine = 1e-6; % m
%%%%%%%%%%%%%%%%%%%%%
x_fine = x_min:spacing_fine:x_max;
fine_npts = length(x_fine);
fine_pts_per_coarse_cell = (fine_npts-1)/(coarse_npts-1) + 1;

% boundary conditions
%%%%%%%%%%%%%%%%%%%%%
T_left = 298;       % Kelvin
T_right = 1000;     % Kelvin
%%%%%%%%%%%%%%%%%%%%%

% define a time step
%%%%%%%%%%%%%%%%%%%%%
dt = 1e-7;      % sec
t_end = 1e-3;      % sec
%%%%%%%%%%%%%%%%%%%%%

% initial conditions
%%%%%%%%%%%%%%%%%%%%%
%T_0 = T_left*ones(coarse_npts,1);   % sudden jump at right side
%T_0(end) = T_right;                 % sudden jump at right side
T_0 = (T_right-T_left)*(0:length(x_coarse)-1)'...   %linear distro
        /length(x_coarse) + T_left;                 %linear distro
%%%%%%%%%%%%%%%%%%%%%


%% generate or load a porosity map on the fine grid
porosity_map_1d = generate_porosity_map_1d(x_fine);
%load('porosity_map_1d_100M.mat');

% generate conductivity from the porosity map by a weighted average
cond_eff_filler = (conductivity_filler/(density_filler*specific_heat_filler));
cond_eff_material = (conductivity_material/(density_material*specific_heat_material));
b = log10(cond_eff_material/cond_eff_filler);
k_fine = cond_eff_filler*10.^(b*porosity_map_1d);
%k_fine = ones(fine_npts,1);


%% Steady state problem
% create the fine matrices Bf, Cf, and qf
% Bf is a tridiagonal: 
tic();
Bf = gallery('tridiag',fine_npts,1,4,1)*spacing_fine/6;
% Cf is also tridiagonal, but the components have contributions from the
% conductivity data
Cf = gallery('tridiag',fine_npts,-1,1,-1)/(2*spacing_fine);
for i=2:fine_npts-1
    Cf(i,i) = Cf(i,i)*(k_fine(i-1)+k_fine(i+1)+2*k_fine(i));
    Cf(i,i-1) = Cf(i,i-1)*(k_fine(i)+k_fine(i-1));
    Cf(i,i+1) = Cf(i,i+1)*(k_fine(i)+k_fine(i+1));
end
% take care of boundaries for Cf
Cf(1,1) = Cf(i,i)*(k_fine(i-1)+k_fine(i+1)+2*k_fine(i));
Cf(1,2) = Cf(i,i+1)*(k_fine(i)+k_fine(i+1));
Cf(fine_npts,fine_npts) = Cf(i,i)*(k_fine(i-1)-k_fine(i+1));
Cf(fine_npts,fine_npts-1) = Cf(i,i+1)*(k_fine(i)+k_fine(i+1));
BC_run_time = toc();

% calculate the D transform matrix by solving subcells 
tic();
D = zeros(coarse_npts, fine_npts);
for i=2:coarse_npts-1;
    f_range_l = ((fine_pts_per_coarse_cell-1)*(i-2)+1):((fine_pts_per_coarse_cell-1)*(i-1)+1);
    f_range_r = ((fine_pts_per_coarse_cell-1)*(i-1)+1):((fine_pts_per_coarse_cell-1)*(i)+1);

    [dl] = fem_steady_homog_diffusion_1d(Cf(f_range_l,f_range_l),...
                                      0, 1);
    [dr] = fem_steady_homog_diffusion_1d(Cf(f_range_r,f_range_r),...
                                      1, 0);
    D(i,f_range_l) = dl';
    D(i,f_range_r) = dr';
    
    
end
% deal with boundaries
f_range_r = 1:(fine_pts_per_coarse_cell+1);
[dr] = fem_steady_homog_diffusion_1d(Cf(f_range_r,f_range_r),...
                                      1, 0);
D(1, f_range_r) = dr';
f_range_l = (fine_npts-fine_pts_per_coarse_cell):fine_npts;
[dl] = fem_steady_homog_diffusion_1d(Cf(f_range_l,f_range_l),...
                                      0, 1);
D(coarse_npts,f_range_l) = dl';
D_run_time = toc();

% calculate the coarse matrices B, C, and q

tic();
B = D*Bf*D';
C = D*Cf*D';
% solve the steady problem now
T_steady = fem_steady_homog_diffusion_1d(C, T_left, T_right);
coarse_run_time = toc();
coarse_run_time = coarse_run_time + D_run_time + BC_run_time;


% solve the fine problem for comparison
tic();
T_fine = fem_steady_homog_diffusion_1d(Cf, T_left, T_right);
fine_run_time = toc();
fine_run_time = fine_run_time + BC_run_time;

%Af = Bf\Cf;
%qf = -Cf(:,1)*T_left -Cf(:,end)*T_right;
%qf = Bf\qf;
%eigmin = eigs(Af, 1, 'sm')
%eigmax = eigs(Af, 1, 'lm')
%stiffness_fine = real(eigmax)/real(eigmin);

%% Output from Steady State problem

%%%%%%%%%%% Plots %%%%%%%%%%%
figure()
hold on
%plot(100*x_fine,porosity_map_1d,'b')
plot(100*x_fine,k_fine,'r');
xlabel('X (cm)')
ylabel('Relative Conductivity m^2/(s-K)')
%legend('porosity map','conductivity map')

figure()
hold on
plot(100*x_fine, T_fine,'r')
plot(100*x_coarse, T_steady,'b')
legend('fine solution','coarse solution')
xlabel('X (cm)')
ylabel('Temp (K)')
title(['Steady Temp. Dist. for (h_c)/(h_f)=',num2str(spacing_coarse/spacing_fine)])


figure()
hold on
whichone = 9;
plot([[0:fine_pts_per_coarse_cell-1]/(fine_pts_per_coarse_cell-1),[(fine_pts_per_coarse_cell-2):-1:0]/(fine_pts_per_coarse_cell-1)],'r')
plot(D(whichone,((whichone-2)*(fine_pts_per_coarse_cell-1)+1):((whichone)*(fine_pts_per_coarse_cell-1)+1)),'k')
title('Modified coarse basis function')
xlabel('X value within basis cell')
ylabel('Basis function value')
legend('Regular hat function','Modified basis function')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Error Calcs and Output %%%%%%%%
disp('%%%%%%%% MsFEM Output %%%%%%%%%%')
disp(['Fine mesh size: ',num2str(spacing_fine)])
disp(['Coarse mesh size: ',num2str(spacing_coarse)])

coarse_pts = 1:fine_pts_per_coarse_cell-1:fine_npts;
error = norm(T_fine(coarse_pts)-T_steady)/norm(T_fine(coarse_pts));
disp(['Error relative to fine mesh is: ',num2str(error)])

disp(['Fine mesh run time: ',num2str(fine_run_time)])
disp(['Coarse mesh run time: ',num2str(coarse_run_time)])
disp('\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Unsteady Problem
% invert the unsteady term matrix 
A = B\C;
% determine the eigenvalues of the problem
eigvals = eig(A);
% evaluate the stiffness
stiffness = max(real(eigvals))/min(real(eigvals));

% deal with boundaries for the reduced problem

q = C(:,1)*T_left +C(:,end)*T_right;
q_red = q(2:end-1);
C_red = C(2:end-1,2:end-1);
B_red = B(2:end-1,2:end-1);

% solve the unsteady problem using implicit scheme
%T_0 = T_steady + 0.05*rand(size(T_steady))*norm(T_steady); % perturbed
%steady state
%
% doing the reduced problem
[times,Y_red] = bwd_euler(B_red, C_red, q_red, dt, 0:dt:t_end, T_0(2:end-1));
Y = [T_left*ones(length(times),1),Y_red,T_right*ones(length(times),1)];
%}
%{
% Using ODE45
odeopts = odeset('Mass',@(t,y) mass_for_ode45(t,y,B));
[times,Y] = ode45(@(t,y) heat_for_ode45(t,y,C,q), 0:dt:t_end, T_0, odeopts);
%}
%[times,Y] = fwd_euler(B, C, q, dt, 0:dt:t_end, T_0); Using forward euler
%[times,Y] = bwd_euler(B, C, q, dt, 0:dt:t_end, T_0); Using backward euler
%[times,Y] = crank_nicolson(B, C, q, dt, 0:dt:t_end, T_0); Using crank
%nicolson



% plot it
figure()
%plot(x_coarse,Y(3,:));
mesh(Y(:,:))
%}



