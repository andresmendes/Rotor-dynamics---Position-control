%% Rotor dynamics - Position control
% Simulation and animation of a controlled rotor with position tracking and
% disturbance rejection.
%
%%

clear ; close all ; clc

%% Circle
th_c    = 0:0.01:2*pi+0.1;      % angle for sweep               [rad]
r       = 1;                    % Radius                        [m]
x_c     = r*cos(th_c);          % Position x                    [m]
y_c     = r*sin(th_c);          % Position y                    [m]

%% Simulation
tf      = 60;                   % Final time                    [s]
fR      = 30;                   % Frame rate                    [fps]
dt      = 1/fR;                 % Time resolution               [s]
time    = linspace(0,tf,tf*fR); % Time                          [s]

% Parameters
% Controller
Kp  = 23;                   % Proportional gain
Ki  = 10;                  % Integral gain
Kd  = 8;                   % Derivative gain
N   = 50;                  % Derivative filter coefficient
sys = @rotor_dynamics;

% Initial conditions
z_c_0   = [0 0];        % Controller dynamics initial conditions
z0      = [z_c_0 0 0];

options = odeset('RelTol',1e-6);
[T, Z]  = ode45(@(t,z) simulation(t,z,sys,Kp,Ki,Kd,N), time, z0,options);

%% Retrieving states
th      = Z(:,3);               % angular position              [rad]
w       = Z(:,4);               % angular speed                 [rad/s]
% Preallocating
torque      = zeros(1,length(T));
disturbance = zeros(1,length(T));
th_ref      = zeros(1,length(T));
for i=1:length(T)
  [dz, torque(i), disturbance(i), th_ref(i)] = simulation(T(i),Z(i,:)',sys,Kp,Ki,Kd,N);
end

% Rotation indicator
x = r*cos(th);                  % Coordinate x                  [m]
y = r*sin(th);                  % Coordinate y                  [m]

%% Animation

figure
set(gcf,'Position',[270   140   640     360  ])

% Create and open video writer object
v = VideoWriter('rotor_dynamics_position.avi');
v.Quality = 100;
open(v);

for i=1:length(time)
    subplot(2,2,1)
      cla
      hold on ; grid on ; axis equal ; box on
      set(gca,'xlim',[-1.2*r 1.2*r],'ylim',[-1.2*r 1.2*r])
      set(gca,'XTick',[], 'YTick', [])
      plot(x_c,y_c,'k','linewidth',3)               % Circle
      plot([-x(i) x(i)],[-y(i) y(i)],'r--')         % Rotation indicator
      plot(0,0,'k*')                                % Origin
      title('Rotor')
    subplot(2,2,2)
      cla
      hold on ; grid on ; box on
      set(gca,'xlim',[T(1) T(end)],'ylim',[min(torque)-0.1*max(torque) max(torque)+0.1*max(torque)])
      plot(T,torque,'b','LineWidth',1.5)
      plot([T(i) T(i)],[min(torque)-0.1*max(torque) max(torque)+0.1*max(torque)],'k--')
      xlabel('Time [s]')
      ylabel('Torque [N.m]')
    subplot(2,2,3)
      cla
      hold on ; grid on ; box on
      set(gca,'xlim',[T(1) T(end)],'ylim',[min(disturbance)-0.1*max(disturbance) max(disturbance)+0.1*max(disturbance)])
      plot(T,disturbance,'b','LineWidth',1.5)
      plot([T(i) T(i)],[min(disturbance)-0.1*max(disturbance) max(disturbance)+0.1*max(disturbance)],'k--')
      xlabel('Time [s]')
      ylabel('Disturbance [N.m]')
    subplot(2,2,4)
      cla
      hold on ; grid on ; box on
      set(gca,'xlim',[T(1) T(end)],'ylim',[min(th*180/pi)-0.1*max(th*180/pi) max(th*180/pi)+0.1*max(th*180/pi)])
      plot(T,th*180/pi,'b','LineWidth',1.5)
      plot(T,th_ref*180/pi,'k')
      plot([T(i) T(i)],[min(th*180/pi)-0.1*max(th*180/pi) max(th*180/pi)+0.1*max(th*180/pi)],'k--')
      xlabel('Time [s]')
      ylabel('Angular position [deg]')
      title('Black=Reference, Blue=Actual')
  
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);

%% Auxiliary functions

function [dz, u, disturbance, R] = simulation(t,z,sys,Kp,Ki,Kd,N)
    % simulation - Simulation of controller and process dynamics.
    %
    % t - time vector
    % 
    % z - States: z1 and z2 are the controller states, z3 is the controlled
    % variable, and so on.
    %
    % sys - Function handle, state space of transfer function.
    %
    % Kp, Ki, Kd - Controller gains
    %
    % N - Derivative filter coefficient
    
    if t < 25
    disturbance = 0;
    else
        disturbance = 10;
    end
    
    % Reference
    if t < 10
        R = pi/2;
    elseif t < 35
        R = 0;
    elseif t < 45
        R = -pi/2;
    else
        R = 0;
    end
    
    e = R - z(3);           % Error

    % Controller derivatives and output
    [dz_c,u]    = pid_controller(z(1:2),e,Kp,Ki,Kd,N);
    % Saturation
    u_max = 20; % Max. torque +- [N.m]
    if u > u_max
        u = u_max;
    end
    if u < -u_max
        u = -u_max;
    end
    % Process derivatives
    dz_sys  = sys(t,z(3:end),u,disturbance);
    % Controller and process derivatives
    dz      = [dz_c ; dz_sys];
end

function dz = rotor_dynamics(~,z,u,disturbance)
    % Parameters
    I   = 2;                    % Moment of inertia             [kg.m2]
    c   = 5;                    % Viscous friction constant     [N.m.s/rad]
    
%     th      = z(1);             % Angular position              [rad]
    w       = z(2);             % Angular speed                 [rad/s]

    % Rotor dynamics
    dz(1,1) = w;
    dz(2,1) = ((u - disturbance) - c*w)/I;
end

function [dz,u] = pid_controller(z,e,Kp,Ki,Kd,N)
    % pid_controller - Controller dynamics and output.
    
    % Constants
    C1 = (Kp/N + Kd);
    C2 = Kp + Ki/N;
    C3 = Ki;

    % State space matrices
    A = [0 1 ; 0 -N];
    B = [0 ; N];
    C = [C3 (C2-C1*N)];
    D = C1*N;

    % State equation
    dz  = A*z + B*e;
    % Output equation
    u   = C*z + D*e;
end