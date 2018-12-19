function liste_out = create_linspace(liste_in, n)
somme = 0;
for i = 1:size(liste_in,1)
    somme = somme + liste_in(i,2)-liste_in(i,1);
end
liste_out=[];
n_total = 0;
for i = 1:size(liste_in,1)
    if i~= size(liste_in,1)
        n_tmp=round(n*(liste_in(i,2)-liste_in(i,1))/somme);
        n_total=n_total+n_tmp;
    else 
        n_tmp=n-n_total;
    end
    liste_out = [liste_out linspace(liste_in(i,1),liste_in(i,2),n_tmp)];
end
end

