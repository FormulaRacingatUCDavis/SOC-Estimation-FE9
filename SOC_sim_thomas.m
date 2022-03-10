%% 
clc; close all; clear
%% Constants
dt = 1;   %Sampling Period
R0 = 0.01;
Rc = 0.015;
Cc = 2400;
Cbat = 18000;
Voc0 = 3.435;
alp = 0.65;

% Noise stuff
Rk = 1e-4;
Qk = 2.5e-7;

%% State Equation
f = @(x, I) [0, 0; 0, -1/(Rc*Cc)]*x + [-1/Cbat; 1/Cc] * I;
h = @(x, I) [alp, -1]*x + -R0*I;

n = 1e3;
x = zeros(2, n);
x(1,1) = 1;
x(2,1) = 0;

xhat = zeros(2, n);
xhat(1,1) = 1;
xhat(2,1) = 0;

y = zeros(1, n);
y(1) = 3.435;
t = zeros(1, n);
dt = 0.01;
uk = 100;
u = zeros(1, n);
Pk = zeros(2,2);

% Continuous-time state-space
Ac = [0, 0; 0, -1/(Rc*Cc)];
Bc = [-1/Cbat; 1/Cc];
Cc = [alp, -1];
Dc = -R0;

% Discretize
aug = [Ac, Bc; zeros(1, 3)];
expmaug = expm(aug * dt);
Ad = expmaug(1:2, 1:2);
Bd = expmaug(1:2, 3);
Cd = Cc;
Dd = Dc;

for i = 1:n 
    if t(i) > 3 && t(i) < 6
        u(1, i) = 10*uk + randn(1, 1) * sqrt(Qk);   
    elseif t(i) > 6 && t(i) < 9
        u(1, i) = -2*uk + randn(1, 1) * sqrt(Qk);
    else
        u(1, i) = uk + randn(1, 1) * sqrt(Qk);
    end
    y(1, i) = h(x(:, i), u(1,i)) + randn(1, 1) * sqrt(Rk);

    if i ~= n
        x(:, i+1) = euler_integration_step(f, x(:, i), u(1, i), dt);
        t(i+1) = t(i) + dt;
        % Kalman filter
        [xhat(:, i+1), Pk] = kalman_filter(Ad, Bd, Cd, Dd, xhat(:, i), u(1, i), y(1, i), Pk, Qk*[1, 0; 0, 0], Rk);
    end
end

subplot(4, 1, 1)
plot(t, x(1, :), 'DisplayName', "SOC")
legend
hold on
plot(t, xhat(1, :), 'DisplayName', "SOChat")
legend
subplot(4, 1, 2)
plot(t, x(2, :), 'DisplayName', "Vc")
hold on
plot(t, xhat(2, :), 'DisplayName', "Vchat")
legend
subplot(4, 1, 3)
plot(t, y(1, :), 'DisplayName', "V")
legend
subplot(4, 1, 4)
plot(t, u(1, :), 'DisplayName', "I")
legend

%% Functions
function [xkp1] = euler_integration_step(f, xk, uk, dt)
    xkp1 = xk + f(xk, uk) * dt;
end

function [xk, Pk] = kalman_filter(Ad, Bd, Cd, Dd, xk, uk, yk, Pk, Qk, Rk)
    % Propagation
    xk = Ad*xk + Bd*uk;
    Pk = Ad*Pk*Ad.' + Qk;

    % Kalman gain
    Lk = Pk*Cd.'*inv(Cd*Pk*Cd.' + Rk);

    % Correction
    xk = xk + Lk*(yk - Cd*xk - Dd*uk);
    Pk = Pk - Lk*Cd*Pk;

end
