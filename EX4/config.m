function config(T)
    for k=1:size(T,2)
        filename=sprintf('config/configuration%d.in', k);
        if exist(filename,'file')~=2
            error(['Le fichier ' filename ' n''existe pas']);
        else
            fileID=fopen(filename);
            C = textscan(fileID,'%s %s','Delimiter',{' ','='},'MultipleDelimsAsOne',1,'EndOfLine','\r\n','CommentStyle','%');
            fclose(fileID);
            fileID=fopen(filename,'w');
            fprintf(fileID,'%% Parametres physiques \n');
            try
                fprintf(fileID,'%-9s = %s \n','nom',T.Properties.VariableNames{k});
                for i=1:size(T,1)
                    if i==10
                        fprintf(fileID,'\n%% Parametres numeriques \n');
                    end
                    if ischar(T{i,k})
                        fprintf(fileID,'%-9s = %s \n',T.Row{i},T{i,k});    
                    else
                        if isnan(T{i,k})
                            disp(C{2}{i});
                            disp(T{i,k});
                            T{i,k}=C{2}{i};
                            disp(T{i,k});
                            fprintf(fileID,'%-9s = %s \n',T.Row{i},T{i,k});
                        else 
                            fprintf(fileID,'%-9s = %.15g \n',T.Row{i},T{i,k});
                        end
                    end
                end
            catch
                warning("erreur lors de l'execution de change_config, rien n'a été changé !");
                fprintf(fileID,'%% Parametres physiques \n');
                for i=1:size(C{1},1)
                    if i==10
                        fprintf(fileID,'\n%% Parametres numeriques \n');
                    end
                    fprintf(fileID,'%-9s = %s \n',C{1}{i},C{2}{i});
                end
            end
            fclose('all');
        end
    end
    disp(T);