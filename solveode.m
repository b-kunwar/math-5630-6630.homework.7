function solveode()
    % Define the ODE as y' = f(t, y)
    f = @(t, y) -2 * t * y; % Example: ODE y' = -2ty
    
    % Initial conditions
    t0 = 0;
    y0 = 1; % y(0) = 1
    h = 0.1; % Step size
    t_end = 2; % Final time
    
    % Time vector
    t = t0:h:t_end;
    n_steps = length(t);
    
    % Initialize solutions
    y_forward_euler = zeros(1, n_steps);
    y_midpoint = zeros(1, n_steps);
    y_rk4 = zeros(1, n_steps);
    
    % Set initial condition
    y_forward_euler(1) = y0;
    y_midpoint(1) = y0;
    y_rk4(1) = y0;
    
    % Forward Euler Method
    for n = 1:n_steps-1
        y_forward_euler(n+1) = y_forward_euler(n) + h * f(t(n), y_forward_euler(n));
    end
    
    % Midpoint Euler Method
    for n = 1:n_steps-1
        k1 = f(t(n), y_midpoint(n));
        y_midpoint(n+1) = y_midpoint(n) + h * f(t(n) + h/2, y_midpoint(n) + (h/2)*k1);
    end
    
    % 4th Order Runge-Kutta Method
    for n = 1:n_steps-1
        k1 = f(t(n), y_rk4(n));
        k2 = f(t(n) + h/2, y_rk4(n) + (h/2)*k1);
        k3 = f(t(n) + h/2, y_rk4(n) + (h/2)*k2);
        k4 = f(t(n) + h, y_rk4(n) + h*k3);
        y_rk4(n+1) = y_rk4(n) + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
    
    % Plot results
    figure;
    plot(t, y_forward_euler, '-o', 'DisplayName', 'Forward Euler');
    hold on;
    plot(t, y_midpoint, '-s', 'DisplayName', 'Midpoint Euler');
    plot(t, y_rk4, '-d', 'DisplayName', 'RK4');
    legend;
    xlabel('t');
    ylabel('y');
    title('Numerical Solutions of ODE');
    grid on;
    hold off;
end
