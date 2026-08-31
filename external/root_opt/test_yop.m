function test_yop(make_plot)
%% Bryson-Denham Example
% Original script authored by Dennis Edblom.
% Adapted from: https://www.yoptimization.com/brysonDenham
%
% Modifications may have been made for integration with this repository.
% All credit for the original formulation and implementation belongs to the author.

bdSystem = YopSystem(...
    'states', 2, ...
    'controls', 1, ...
    'model', @trolleyModel ...
    );

time = bdSystem.t;
trolley = bdSystem.y;

ocp = YopOcp();
ocp.min({ timeIntegral( 1/2*trolley.acceleration^2 ) });
ocp.st(...
    'systems', bdSystem, ...
    ... % Initial conditions
    {  0  '==' t_0( trolley.position ) }, ...
    {  1  '==' t_0( trolley.speed    ) }, ...
    ... % Terminal conditions
    {  1  '==' t_f( time ) }, ...
    {  0  '==' t_f( trolley.position ) }, ...
    { -1  '==' t_f( trolley.speed    ) }, ...
    ... % Constraints
    { 1/9 '>=' trolley.position        } ...
    );

% Solving the OCP
sol = ocp.solve('controlIntervals', 20, ...
    'ipopt', struct('max_iter', 2000, "print_level",1, "max_cpu_time", 500,...
    'tol',1e-6));

% Plot the results
if(make_plot==true)
    figure(1)
    subplot(211); hold on
    sol.plot(time, trolley.position)
    xlabel('Time')
    ylabel('Position')

    subplot(212); hold on
    sol.plot(time, trolley.speed)
    xlabel('Time')
    ylabel('Velocity')

    figure(2); hold on
    sol.stairs(time, trolley.acceleration)
    xlabel('Time')
    ylabel('Acceleration (Control)')
end
if(sol.converged==1)
 % If no error is thrown, the test passes
    disp('test_yop passed.');
end
end
%%
function [dx, y] = trolleyModel(time, state, control)

position = state(1);
speed = state(2);
acceleration = control;
dx = [speed; acceleration];

y.position = position;
y.speed = speed;
y.acceleration = acceleration;

end
