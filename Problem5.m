Nvals = [20 40 80 160];
L1errors = [0 0 0 0];
LInferrors = [0 0 0 0];
L2errors = [0 0 0 0];
H1errors = [0 0 0 0];

for j = 1:4
    N = Nvals(j);
    h = 1/(N+1);
    
    f = @(x, y) 2*pi*pi*sin(pi*x).*cos(pi*y);
    uexact = @(x, y) sin(pi*x).*cos(pi*y);
    
    x = h*(0:N+1);
    y = h*(N+1:-1:0);
    [plotX, plotY] = meshgrid(x, y);
    
    xvals = h*(1:N);
    yvals = h*(N:-1:1);
    [X, Y] = meshgrid(xvals, yvals);
    
    tempF = h*h*f(X, Y);
    tempF(1,:) = tempF(1,:) - sin(pi*xvals);
    tempF(N,:) = tempF(N,:) + sin(pi*xvals);
    
    F = reshape(fliplr(tempF.'), N*N, 1);
    
    maindiag = 4*ones(N*N, 1);
    subdiag1 = -1*ones((N*N)-1, 1);
    for k = 1:N-1
        subdiag1(k*N) = subdiag1(k*N) + 1;
    end
    subdiag2 = -1*ones(N*(N-1), 1);
    
    A = diag(maindiag) + diag(subdiag1, -1) + diag(subdiag1, 1) + diag(subdiag2, N) + diag(subdiag2, -N);
    
    Uvec = A \ F;
    
    U = fliplr(reshape(Uvec, N, N)).';
    U = [zeros(N, 1) U zeros(N, 1)];
    U = [-sin(pi*x); U; sin(pi*x);];
    
    Utrue = uexact(plotX, plotY);
    
    error = Utrue - U;
    
    L1errors(j) = max(sum(abs(error)));
    LInferrors(j) = norm(error, "inf");
    L2errors(j) = norm(error, 2);

    dxErrors = [diff(error, 1, 2) / h, (error(:, end) - error(:, end-1)) / h];
    dyErrors = [diff(error, 1, 1) / h; (error(end, :) - error(end-1, :)) / h];
    gradErrors = dxErrors.^2 + dyErrors.^2;

    H1errors(j) = L2errors(j) + norm(gradErrors, 2);

end

L1errors
L2errors
LInferrors
H1errors

hvals = (Nvals + 1).^-1

coeff = polyfit(log(hvals), log(L2errors), 1);
coeff(1)

loglog(hvals, L1errors, "LineWidth", 2)
hold on;
loglog(hvals, L2errors, "LineWidth", 2)
hold on;
loglog(hvals, LInferrors, "LineWidth", 2)
hold on;
loglog(hvals, H1errors, "LineWidth", 2)
xlabel("log of the stepsize")
ylabel("log of the error")
title("loglog Plot of Error vs Stepsize for the Elliptic BVP")
legend("L^{1} Norm", "L^{2} Norm", "L^{inf} Norm", "H^{1} Norm")

%{
surf(plotX, plotY, U)
%surf(plotX, plotY, uexact(plotX, plotY))
xlabel("x")
ylabel("y")
zlabel("f(x,y)")
%}