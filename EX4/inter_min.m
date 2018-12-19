function [v,i]=inter_min(x,y,n)
    [value,indice]=min(y);
    x=x(indice-n:indice+n);
    y=y(indice-n:indice+n);
    [p,~,mu]=polyfit(x,y,n);
    x=linspace(x(1),x(end),100000);
    y=polyval(p,x,[],mu);
    [v,i]=min(y);
    i=x(i);
end
