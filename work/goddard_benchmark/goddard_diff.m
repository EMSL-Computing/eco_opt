%% Comparison between Yop, CasADi, PROPT
load('yop_goddard.mat'); 
load('casadi_goddard.mat');
load('propt_goddard.mat');

figure(1)
subplot(311); hold on
plot(yop_t,    yop_x(:,1), 'r', 'LineWidth', 2)
plot(casadi_t, casadi_x(:,1), 'k')
plot(propt_t,  propt_x(:,1), 'b')
legend('Yop', 'CasADi', 'PROPT')
xlabel('Time'); ylabel('Velocity')

subplot(312); hold on
plot(yop_t,    yop_x(:,2), 'r', 'LineWidth', 2)
plot(casadi_t, casadi_x(:,2), 'k')
plot(propt_t,  propt_x(:,2), 'b')
legend('Yop', 'CasADi', 'PROPT')
xlabel('Time'); ylabel('Height')

subplot(313); hold on
plot(yop_t,    yop_x(:,3), 'r', 'LineWidth', 2)
plot(casadi_t, casadi_x(:,3), 'k')
plot(propt_t,  propt_x(:,3), 'b')
legend('Yop', 'CasADi', 'PROPT')
xlabel('Time'); ylabel('Mass')

figure(2); hold on
stairs(yop_t,    yop_u, 'r', 'LineWidth', 2)
stairs(casadi_t, casadi_u, 'k')
stairs(propt_t,  propt_u, 'b')
legend('Yop', 'CasADi', 'PROPT')
xlabel('Time'); ylabel('F (Control)')