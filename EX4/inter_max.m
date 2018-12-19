function [v,i]=inter_max(x,y,n)
    [v,i]=inter_min(x,-y,n);
    v=-v;
end