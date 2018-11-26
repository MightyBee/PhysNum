function C=change_parameter(configfileNb,parameterName,parameterValue)

    filename=sprintf('config/configuration%d.in', configfileNb);
    if exist(filename,'file')~=2
        error(['Le fichier ' filename ' n''existe pas']);
    else
        fileID=fopen(filename);
        C = textscan(fileID,'%s %s','Delimiter',{' ','='},'MultipleDelimsAsOne',1,'EndOfLine','\r\n','CommentStyle','%');
        fclose(fileID);
        fileID=fopen(filename,'w');
        fprintf(fileID,'%% Parametres physiques \n');
        try
            for i=1:size(C{1},1)
                if ischar(parameterValue) && strcmp(parameterName,C{1}{i})
                    C{2}{i}=parameterValue;
                elseif strcmp(parameterName,C{1}{i})
                    C{2}{i}=sprintf('%.15g', parameterValue);
                end
                if (configfileNb==0 && i==6) || (configfileNb>0 && i==9)
                    fprintf(fileID,'\n%% Parametres numeriques \n');
                end
                fprintf(fileID,'%-9s = %s \n',C{1}{i},C{2}{i});
            end
        catch
            warning("erreur lors de l'execution de change_config, rien n'a été changé !");
            fprintf(fileID,'%% Parametres physiques \n');
            for i=1:size(C{1},1)
                if i== 10
                    fprintf(fileID,'\n%% Parametres numeriques \n');
                end
                fprintf(fileID,'%-9s = %s \n',C{1}{i},C{2}{i});
            end
        end
        fclose('all');
    end
