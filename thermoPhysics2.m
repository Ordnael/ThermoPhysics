clear; close all; clc;

tnow = datestr(now,30);
N = 10000;
E_limit = 1e-12;

% Vectors.
dx = 0.0001;
x = [0:dx:0.002];
y = [0:dx:0.001];
y = fliplr(y);  % To start from top.
y = y';         % y as a column vector.

% Boundary conditions.
Ta = 20;
hmin = 2000;
hmax = 8000;
hx = linspace(hmax,hmin,size(x,2));
k = 16;
qv = 2e8;

% Temperature matrix.
s1 = size(x,2);
s2 = size(y,1);
T = zeros(size(y,1)+2, size(x,2)+2);
E = zeros(N,1);
T = T + 100;    % Initial T distribution.

% Applying BC at the borders.

% Creating the buffer matrix.
TT = T;
it = 1;         % Iteration counter.
E_act = 1;

while E_act > E_limit
    % Break if iteration limit reached.
    if it > N
        STOP = 'Iteration limit reached!'
        break
    end
    
    % Loop over the horizontal external nodes.
    for n = 2:s1+1
        m = 1;
        % Upper boundary.
        TT(m,n) = 2*dx*hx(1,n-1)/k*(Ta-T(m+1,n)) + T(m+2,n);
        % Lower boundary.
        m = s2+2;
        TT(m,n) = T(m-2,n);
    end
    
    % Loop over the vertical external nodes.
    for m = 2:s2+1
        n = 1;
        TT(m,n) = T(m,n+2);
        n = s1+2;
        TT(m,n) = 2*dx*hmin/k*(Ta-T(m,n-1))+T(m,n-2);
    end
    
    
    % Loop over the internal nodes.
    for n=2:s1+1
        for m=2:s2+1
            TT(m,n) = 0.25*(T(m+1,n)+T(m-1,n)+T(m,n+1)+T(m,n-1)) + qv*dx^2/4/k;
        end
    end
    
    % Error calculation.
    E(it) = max(max(abs(T(2:s2+1,2:s1+1)-TT(2:s2+1,2:s1+1))));
    E_act = E(it);
    
    % Prepare for next iteration.
    T = TT;
    it = it + 1;
end

[X,Y] = meshgrid(x,y);
h1 = figure(1);
contour(X,Y,T(2:s2+1,2:s1+1),40);
xlabel('x, m');
ylabel('y, m');
zlabel('T, °C');
colorbar;
grid on;
saveas(h1,[tnow '_temp_contour.fig']);
print('-dpng',[tnow '_temp_contour']);

h2 = figure(2);
semilogy(E);
xlabel('Number of iterations, -');
ylabel('Maximum error, °C');
grid on;
saveas(h2,[tnow '_error.fig']);
print('-dpng',[tnow '_error']);

h3 = figure(3);
plot(x, T(round(size(T,1)/2),2:s1+1));
legend(['y = ' num2str(y(round(size(T,1)/2),1)) 'm']);
xlabel('x, m');
ylabel('y, °C');
saveas(h3,[tnow '_Tx.fig']);
print('-dpng',[tnow '_Tx']);

p4.IN.x = x;
p4.IN.y = y;
p4.IN.dx = dx;
p4.IN.Ta = Ta;
p4.IN.hmin = hmin;
p4.IN.hmax = hmax;
p4.IN.hx = hx;
p4.IN.k = k;
p4.IN.qv = qv;
p4.IN.N = N;
p4.IN.E_limit = E_limit;

p4.OUT.T = T;
p4.OUT.E = E;
p4.OUT.it = it;

save([tnow '_p4'], 'p4');