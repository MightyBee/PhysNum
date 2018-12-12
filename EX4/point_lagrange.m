function x0=point_lagrange(xM,xP,epsilon)
    yM=f(xM);
    yP=f(xP);
    fun = @(x) f(x);
    disp(sprintf('%.15g',fzero(fun,[xM xP])));
    x0=xM-yM*(xP-xM)/(yP-yM);
    y0=f(x0);
    while abs(y0)>epsilon
        if sign(yP-yM)*y0 > 0
            xP=x0;
            yP=y0;
        else
            xM=x0;
            yM=y0;
        end
        x0=xM-yM*(xP-xM)/(yP-yM);
        y0=f(x0);
    end


function y=f(x)
    d=384748000;
    mT=5.972e24;
    mL=7.3477e22;
    alpha=mL/(mT+mL);
    beta=mT/(mT+mL);
    rG=alpha*d;
    xT=0-rG;
    xL=d-rG;
    y=d^3*(beta*(x-xT)./abs(x-xT).^3+alpha*(x-xL)./abs(x-xL).^3)-x;
