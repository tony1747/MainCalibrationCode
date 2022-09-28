function outStruct = paramInfo(name, givenInfo)
% PARAMINFO
% This function returns default parameter info for 
%



bool_givenInfo = false;
if exist('givenInfo','var')
    if string(class(givenInfo)) == "struct"
        bool_givenInfo = true; 
    else
        msg = strcat("Error. ", "paramInfo input ''",inputname(1), ...
            "'' of class ''", class(givenInfo), ...
            "'',",newline, " expected class ''struct''.");
        error(msg);
    end
end


paramfile = strcat("paramInfo_",string(name),".csv");
if(isfolder("models") & ~isfile(paramfile))
    paramfile = strcat('models/',paramfile);
end

if(~isfile(paramfile))
    
    msg = strcat("Error: cannot find file: ",paramfile);
    error(msg)
end


filetable = readtable(paramfile, 'ReadRowNames',true);
rowsToFind = ["lb","ub","default","inds"];
if(~all(ismember(rowsToFind,filetable.Row)))
    msg = strcat("Error. Table ''", paramfile, "'' has improper format.");
    error(msg);
end


%%% this check is obsolete
if(ismember("names",filetable.Row))
    names = string(table2array(filetable("names",:)));
else
    names = string(filetable.Properties.VariableNames);
end

N_p = length(names);
outStruct.names = names;
%%% old form, used to use a:b instead of logicals
%outStruct.inds = eval(strjoin(string(filetable{"inds",[1 2]}),":"));
outStruct.inds = find(filetable{"inds",:});

for row=rowsToFind(1:3)
    outStruct.(row) = filetable{row,:};
end
    
%%% add any given information
if(bool_givenInfo)
    for key=intersect(fieldnames(outStruct),fieldnames(givenInfo))
        outStruct.(string(key)) = givenInfo.(string(key));
    end
end


end
