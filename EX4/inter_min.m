function v=inter_min(x,y,n)
    [value,indice]=min(y);
    x=x(indice-n:indice+n);
    y=y(indice-n:indice+n);
    [p,~,mu]=polyfit(x,y,n);
    x=linspace(x(1),x(end),100000);
    y=polyval(p,x,[],mu);
    v=min(y);
end
%% sauvergarde %%

% elseif strcmp(paraName,'condIn')   
%     paramstr = {"vx0"; "vy0"};
%     nsimul=round(sqrt(nsimul));
%     if nsimul<3
%         nsimul=3;
%     end
%     theta0 = linspace(0,2*pi,nsimul+1);
%     theta0 = theta0(1:nsimul);
%     v0A = linspace(10,10000,nsimul);
%     v0A = [v0A(1)*ones(1,nsimul) v0A(2:nsimul-1) v0A(nsimul)*ones(1,nsimul)];
%     theta=[theta0 zeros(1,nsimul*(nsimul-2)) theta0];
%     for i=2:nsimul-1
%        theta((i-1)*nsimul+1:i*nsimul) = theta0;
%        v0A=[v0A(1:2*nsimul-1-i) v0A(2*nsimul-i)*ones(1,nsimul) v0A(2*nsimul-i+1:end)];
%     end
%     param = [v0A.*cos(theta); v0A.*sin(theta)+omega*rA]; % Valeurs du parametre a scanner
%     nsimul = nsimul^2;
%     configfileNb=3;