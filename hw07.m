% Author: Bikash Kunwar / bzk0067@auburn.edu
% Date: 2024-09-01
% Assignment Name: hw07

classdef hw07

    methods (Static)
        function y = p1(func, y0, tspan, n_steps, method)
            % Solves the ODE y' = f(t, y) with initial condition y(t0) = y0 using the specified method (euler, rk4, midpoint).
            % over the interval tspan=[a, b] with n_steps. The function f(t, y) is provided as a function handle.
            %
            %:param func: function handle f(t, y) that defines the ODE y' = f(t, y)
            %:param y0: initial condition y(t0) = y0
            %:param tspan: interval [a, b] over which to solve the ODE
            %:param n_steps: number of steps to take to solve the ODE, interval size = (b-a)/n_steps.
            %:param method: string that specifies the method to use. It can be 'euler', 'midpoint', or 'rk4'
            %
            %:return: none, but plots the solution y(t) over the interval tspan

            % Your implementation here. Euler method is implemented for an example. Implement the other methods.

            t0 = tspan(1);
            tf = tspan(2);

            h = (tf - t0) / n_steps;
            t = t0:h:tf;
            y = zeros(1, length(t));

            y(1) = y0;

            if strcmp(method, 'euler')
                for i = 1:n_steps
                    k1 = func(t(i), y(i));
                    y(i+1) = y(i) + h * k1;
                end
            elseif strcmp(method, 'rk4')
                % your code for Runge-Kutta 4 method here
                for i = 1:n_steps
                    k1 = func(t(i), y(i));
                    k2 = func(t(i) + h/2, y(i) + h/2 * k1);
                    k3 = func(t(i) + h/2, y(i) + h/2 * k2);
                    k4 = func(t(i) + h, y(i) + h * k3);
                    y(i+1) = y(i) + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
                end

            elseif strcmp(method, 'midpoint')
                % your code for Midpoint method here
                for i = 1:n_steps
                    k1 = func(t(i), y(i));
                    y_mid = y(i) + (h/2) * k1; % Midpoint
                    t_mid = t(i) + h/2; % Time at midpoint
                    y(i+1) = y(i) + h * func(t_mid, y_mid); % Full step
                end

            else
                error('Invalid method. Choose "euler", "rk4", or "midpoint".');
            end
        end

        function p2(method)
            % Test the implemented methods on the ODE
            % y' = t(y - t sin(t)) with initial condition y(0) = 1 over the interval [0, 1], [0, 3], and [0, 5] with variable step lengths.

            %
            % Plot the solution y(t) over the interval for various step sizes h. And plot the
            % exact solution y(t) = t sin(t) + cos(t) over the same interval.
            %
            % Run the commands below to test the implemented methods:
            %
            %> hw07.p2('euler');
            %> hw07.p2('midpoint');
            %> hw07.p2('rk4');
            %
            % Observe the solution and error plots for the numerical solutions with different step sizes. Write your observations in the comments.

            % Your comment here (e.g, how does the error change with step size and the time span, etc.):

            % Euler method
            % As step size is decreased the numerical solution approaches the exact solution, however doesn't
            % exactly match the exact solution. As time span is increased,
            % the solution seems to have better agreement for interval
            % [0,3] comared to tspan [0,1], but deviations increase for [0,5]. The convergence is first
            % order.

            % Midpoint Method
            % The oder of convergence is close to second order for all time
            % spans. The agreement between numerical solution and exact solution is  better compared to euler method for all step size for tspan [0,1] and [0,3].
            % However, the disagreement is huge for tspan [0,5], especially for bigger step size, which seem to work fine on tspan [0,1] and [0,3].

            % RK4
            % On the interval [0,1], and [0,3], all step size results close to exact solution. The rate of convergence is 4th order. As step size is decreased, error reduces at that order.
            % On the interval [0,5], the error is reduced greatly close to
            % t = 5, where euler and midpoint method gave larger error. The
            % error magnitude itself dropped to 10^-2 (where we got postive
            % error with the previous methods). The error is O(h^4) for
            % some tspan [0,1] not for [0,3] and [0,5]





            f = @(t, y) t * (y - t * sin(t));

            tf_values = [1,3,5];

            figure('Position', [0 0 1200 2000]);
            for tf_index = 1:length(tf_values)

                t0 = 0; tf = tf_values(tf_index); y0 = 1;
                exact_sol = @(t) t .* sin(t) + cos(t);
                error = zeros(1, 8);
                h = [1e-1, 1e-1 * 10/11, 1e-1*4/5, 1e-1*0.72, 1e-1*0.64, 1e-1*1/2, 1e-1*2/5, 1e-1*0.32];
                c = {'g--', 'g-.', 'b--', 'b-.', 'm--', 'm-.', 'k--', 'k-.'};

                subplot(length(tf_values),2, tf_index * 2 - 1);
                title(['Numerical', ' solution using ', method, ' method and Exact Solutions', ' on interval [', num2str(t0), ', ', num2str(tf), ']']);
                for i = 1:length(h)
                    n_steps = (tf - t0) / h(i);
                    y=hw07.p1(f, y0, [t0, tf], n_steps, method);
                    error(i) = max(abs(y - exact_sol(t0:h(i):tf)));
                    hold on;
                    plot(t0:h(i):tf, y, sprintf('%s', c{i}), 'DisplayName', ['h = ', num2str(h(i))]);
                end

                plot(t0:h(end):tf, exact_sol(t0:h(end):tf), 'r-', 'DisplayName', 'Exact Solution');
                hold off; legend("Location", 'best'); grid on; xlabel('t'); ylabel('y(t)')


                subplot(length(tf_values),2, tf_index * 2);
                loglog(h, error, 'b-o', 'DisplayName', 'Max Error vs. Step Size');
                hold on;
                loglog(h, h.^4 * error(1)/h(1)^4, 'r--', 'DisplayName', '4th order convergence');
                loglog(h, h.^3 * error(1)/h(1)^3, 'g--', 'DisplayName', '3rd order convergence');
                loglog(h, h.^2 * error(1)/h(1)^2, 'm--', 'DisplayName', '2nd order convergence');
                loglog(h, h.^1 * error(1)/h(1)^1, 'k--', 'DisplayName', '1st order convergence');
                hold off;
                title(['Error vs. Step Size for y'' = t(y - t sin(t))', ' on [', num2str(t0), ', ', num2str(tf), ']']); xlabel('Step Size (h)'); ylabel('Error'); grid on; legend("Location", 'best');
            end
        end

        function p3()
            % For 6630 ONLY
            % First implement the 3/8 rule for Runge Kutta method.
            %
            % The implementation should be done in the function rk4_38_rule below. It is a subfunction which can only be called within p3 method.
            % Then run hw07.p3() and compare the results with the 4th order Runge Kutta method. Write your observations in the comments.
            %
            % Your comment here (e.g, how does the error change with step size and the time span, is there a clear difference in the running time and error (you may need to run a few times to conclude), etc.):
            % The 3/8 rule performs slightly better than Rk4 in terms of
            % error. The convergence is 4 th order for 3/8
            % and slightly less for RK4. However, the difference is minimal to make any
            % practical difference. The run times are also similar with minor fluctations. 
            % At very small step sizze (< 10^-3), the
            % error, convergecne and run time are almost coincident.
        

            function y = rk4_38_rule(func, y0, tspan, n_steps)
                % rk4_38_rule: Runge-Kutta method with 4th order and 3/8 rule for a system of ODEs.
                %:param func: function handle f(t, y) that defines the ODE y' = f(t, y)
                %:param y0: initial condition y(t0) = y0
                %:param tspan: interval [a, b] over which to solve the ODE
                %:param n_steps: number of steps to take to solve the ODE, interval size = (b-a)/n_steps.

                t0 = tspan(1);
                tf = tspan(2);

                h = (tf - t0) / n_steps;
                t = t0:h:tf;
                y = zeros(1, length(t));

                y(1) = y0;

                % write your code here.
                % Loop through each time step
         
                for i = 1:n_steps
                    tn = t(i);
                    yn = y(i);

                    % Compute the Runge-Kutta coefficients
                    k1 = func(tn, yn);
                    k2 = func(tn + h/3, yn + h*k1/3);
                    k3 = func(tn + 2*h/3, yn - h*k1/3 + h*k2);
                    k4 = func(tn + h, yn + h*k1 - h*k2 + h*k3);

                    % Update the solution using the 3/8 rule
                    y(i+1) = yn + h/8 * (k1 + 3*k2 + 3*k3 + k4);
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Do not modify the code below this line.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            f = @(t, y) t * (y - t * sin(t));

            tf_values = [3, 5, 7];
            figure('Position', [0 0 1200 2000]);
            for tf_index = 1:length(tf_values)
                t0 = 0; tf = tf_values(tf_index); y0 = 1;
                exact_sol = @(t) t .* sin(t) + cos(t);
                error1 = zeros(1, 8);
                error2 = zeros(1, 8);
                runtime1 = zeros(1, 8);
                runtime2 = zeros(1, 8);
                hs = 2.^(-1:-1:-8) * 1e-1;

                for ind = 1:length(hs)
                    n_steps = (tf - t0) / hs(ind);

                    tic;% run 100 times and take the average time
                    for run = 1:100
                        y1=rk4_38_rule(f, y0, [t0, tf], n_steps);
                    end
                    runtime1(ind) = toc/100;

                    tic;% run 100 times and take the average time
                    for run = 1:100
                        y2=hw07.p1(f, y0, [t0, tf], n_steps, 'rk4');
                    end
                    runtime2(ind) = toc/100;

                    error1(ind) = max(abs(y1 - exact_sol(t0:hs(ind):tf)));
                    error2(ind) = max(abs(y2 - exact_sol(t0:hs(ind):tf)));
                end

                subplot(length(tf_values),2, tf_index * 2 - 1);
                loglog(hs, error1, 'b-o', 'DisplayName', 'Max Error (3/8 Rule) vs. Step Size');
                hold on;
                loglog(hs, error2, 'g-d', 'DisplayName', 'Max Error (RK4) vs. Step Size');
                loglog(hs, hs.^4 * error1(1)/hs(1)^4 , 'r--', 'DisplayName', '4th order convergence');
                hold off;
                title(['Error vs. Step Size for y'' = t(y - t sin(t))', ' on [', num2str(t0), ', ', num2str(tf), ']']); xlabel('Step Size (h)'); ylabel('Error'); grid on; legend("Location", 'best');
                subplot(length(tf_values),2, tf_index * 2);
                loglog(hs, runtime1, 'r-o', 'DisplayName', 'Runtime (3/8 Rule) vs. Step Size');
                hold on;
                loglog(hs, runtime2, 'm-d', 'DisplayName', 'Runtime (RK4) vs. Step Size');
                hold off;
                title(['Runtime vs. Step Size for y'' = t(y - t sin(t))', ' on [', num2str(t0), ', ', num2str(tf), ']']); xlabel('Step Size (h)'); ylabel('Runtime (s)'); grid on; legend("Location", 'best');
            end

        end
    end
end