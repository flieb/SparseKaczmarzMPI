function [l1,l2,c] = refineLambda(lambda1, lambda2, ffun, pval,x)

    cold = pval;
    l2 = lambda2;
    l1 = lambda1;
    iter1 = 0;
    fac = 0.4;
    minfac = 0.01;
    
    while (1)
        
        %optimize psnr by adjusting lambda1 until maximum:
        ffun2 = @(x) ffun(x,l2);
        [l1,iter2,c1] = optim_var(l1,ffun2,fac,cold,x);

        %optimize psnr by adjusting lambda2 until maximum:
        ffun2 = @(x) ffun(l1,x);
        [l2,iter3,c] = optim_var(l2,ffun2,fac,c1,x);

        
        
        iter1 = iter1 + 1;
        fprintf('iter = %d, iter2=%d, iter3=%d, fac=%4.3f, psnr=%4.3f\n',iter1,iter2,iter3,fac,c);
        
        %decrease factor by which lambda is varied in each iteration
        if (iter2==0 && iter3==0 && fac~=minfac)
            fac = max(fac / 2,minfac);
            continue;
        end
        
        if (abs(c-cold)<0.001)
            break;
        end
        cold = c;
        fac = max(fac / 2,minfac);
        
    end
end


function [ll,iter,c] = optim_var(val,ffun,factor,pval,x) 
    iter = 0;
    d = val*factor;
    c_old = pval;
    
    %first check if adding increases pnsr
    l = val + d;
    c = psnr(ffun(l),x);
    while(c>c_old)
        c_old = c;
        l = l + d;
        c = psnr(ffun(l),x);
        iter = iter + 1;
    end
    if (iter>0)
        ll = l - d;
        c = psnr(ffun(ll),x);
        return;
    end
    
    %then subtract 
    l = val - d;
    c = psnr(ffun(l),x);
    while (c>c_old)
        c_old = c;
        l = l - d;
        if (l<0) l = l + d; end
        c = psnr(ffun(l),x);
        iter = iter + 1;
    end
    ll = l + d;
    c = psnr(ffun(ll),x);
end