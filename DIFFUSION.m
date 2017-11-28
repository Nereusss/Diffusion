% NUMERICALLY SIMULATE DIFFUSION EQUATION

% Prepare the Editor, Workspace and Command Window
clc
clear
close all

% SETTING CONSTANTS

k = 1;      % DIFFUSIVITY
dt = 0.00005; % TIME STEP
dx = 0.01;   % SPACIAL STEP(IN THE CASE OF ONE DIMENSION)
c = (k*dt)/(dx*dx); % CFL CONDITION

X = 1;    % UNIT LENGTH OF SPACIAL DOMAIN
T = 0.4;    % WHOLE TIME FOR APPROXIMATION
N = T/dt; % NUMBERS OF TIME STEPS
M = X/dx; % NUMBERS OF SPACIAL STEPS

U = zeros(M,N); % INITIALIZE ALL ZEROES FOR MATRIX

% INITIALIZING CONDITIONS

for k = 2:1:M-1
    U(k,1) = 1+ cos(pi.*(k/M)); % GIVEN BY f(x)=1+cos(pi.*x)
end

% APPROXIMATION

for j = 1:1:N-1
    for i = 2:1:M-1
        U(1,j) = U(2,j); % ZERO BOUNDARY CONDITON
        U(M,j) = U(M-1,j); % ZERO BOUNDARY CONDITION
        U(i,j+1) = U(i,j)+c.*(U(i-1,j)-2.*U(i,j)+U(i+1,j)); % DISCRETE APPROXIMATION OF THE PDE
    end
end

% COMPLETING BOUNDARY CONDITIONS AT THE FINAL TIME STEP 
U(1,N) = U(2,N);
U(M,N) = U(M-1,N);

% PLOTTING THE GRAPH
imagesc(U)
xlabel('Time Steps');
ylabel('The bar from top to end (divided into 100 pieces)');
colorbar