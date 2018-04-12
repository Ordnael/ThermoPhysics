clear; close all; clc;

tnow = datestr(now,30);
N = 10000;
E_limit = 1e-12;

% Vectors.
x = [0:0.001:0.02];
y = [0:0.001:0.01];
y = fliplr(y);  % To start from top.
y = y';         % y as a column vector.

% Initial temperatures.
BCT1 = 0;
BCT2 = 150;
BCT3 = 100;
BCT4 = 100;

% Temperature matrix.
T = zeros(size(y,1), size(x,2));
E = zeros(N,1);
T = T + 100;    % Initial T distribution.

% Applying BC at the borders.
T(1,:) = linspace(BCT1,BCT2,size(x,2));
T(end,:) = linspace(BCT3,BCT4,size(x,2));
T(:,1) = linspace(BCT1,BCT3,size(y,1));
T(:,end) = linspace(BCT2,BCT4,size(y,1));

% Creating the buffer matrix.
TT = T;
it = 1;         % Iteration counter.
E_act = 1;

while E_act > E_limit
    % Break if iteration limit reached.
    if it > N
        STOP = 'Iteration limit reached!';
        break
    end
    
    % Loop over the internal nodes.
    for n=2:size(T,2)-1
        for m=2:size(T,1)-1
            TT(m,n) = 0.25*(T(m+1,n)+T(m-1,n)+T(m,n+1)+T(m,n-1));
        end
    end
    
    % Error calculation.
    E(it) = max(max(abs(T-TT)));
    E_act = E(it);
    
    % Prepare for next iteration.
    T = TT;
    it = it + 1;
end

[X,Y] = meshgrid(x,y);
h1 = figure(1);
contour(X,Y,T,40);
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
plot(x, T(round(size(T,1)/2),:));
legend(['y = ' num2str(y(round(size(T,1)/2),1)) 'm']);
xlabel('x, m');
ylabel('y, °C');
saveas(h3,[tnow '_Tx.fig']);
print('-dpng',[tnow '_Tx']);

p3.IN.x = x;
p3.IN.y = y;
p3.IN.BCT1 = BCT1;
p3.IN.BCT2 = BCT2;
p3.IN.BCT3 = BCT3;
p3.IN.BCT4 = BCT4;
p3.IN.N = N;
p3.IN.E_limit = E_limit;

p3.OUT.T = T;
p3.OUT.E = E;
p3.OUT.it = it;

save([tnow '_p3'], 'p3');