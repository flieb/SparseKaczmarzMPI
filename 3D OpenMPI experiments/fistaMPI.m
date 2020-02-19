function [c,convergence,kk] = fistaMPI(A,b,it,lambda,dim,epsilon,method,gamma,Ap,AtA,Atb)
% implements fista as iterative thresholding
% Copyright (C) 2020 Florian Lieb


plotting = 0;
fignum = 4;

[n,m] = size(A);
u = complex(zeros(m,1));



tn = 1;

if nargin < 10
    gradfun = @(x) x - gamma*Ap*(A*x-b);
else
    gradfun = @(x) x - gamma*(AtA*x - Atb);
end




convergence = zeros(1,it);

for kk=1:it
    y = gradfun(u);
    %tmp = A*u-b;
    %y = u - gamma*Ap*tmp;
    
   
    %choose proximal operator
    y = myProxL1_swt3(real(y),lambda,dim,method);
    %y = myProxL1_PySwtn(real(y),lambda,dim,method);
    
    y(real(y) < 0) = 0;
    
    tn1 = (1 + sqrt(1 + 4*tn^2))/2;
    zs=y+(tn-1)/tn1*(y-u);
    %zs = y;
    
    convergence(kk) = norm(u-zs,'fro')./(norm(u,'fro')+1e-3);
    %fprintf('%5.6f\n',convergence(kk));
    
    %updates
    tn = tn1;
    u = zs;
    
    
    %convergence(kk) = norm(u-zs,'fro')./(norm(u,'fro')+1e-3);
    
    if kk>2
        %if (abs(cngPS(kk)-cngPS(kk-1)) < epsilon)
        if (convergence(kk) < epsilon)
            break
        end
    end
    
    if plotting
        plotc(real(zs),dim,fignum);
        drawnow;
    end

end

%figure, semilogy(convergence);
c = u;