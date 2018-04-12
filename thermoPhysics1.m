clear; close all; clc;

tnow = datestr(now,30);
N = 10000;

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

while E_act > 1e-12
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

surf(T);