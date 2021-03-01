% This function returns the ID given as inputs the path to the
% file containing the output of the ID and the joints of
% interest
%
% Author: Antoine Falisse
% Date: 12/19/2018
%  
% Modified: Aaron Fox
% Date: 09/11/2020
%
function ID = getID(pathID,fileID,joints)

    import org.opensim.modeling.*
    
    %Extract joint kinematics
    %Load id .sto file
    idData = TimeSeriesTable([pathID,'\',fileID]);
    %Get time
    for tt = 0:idData.getIndependentColumn().size()-1
        ID.time(tt+1,1) = idData.getIndependentColumn().get(tt);
    end
    clear tt
    %Get dynamics data
    for jj = 1:length(joints)
        %Check for translational coordinate
        if strcmp(joints{jj},'pelvis_tx') || ...
                strcmp(joints{jj},'pelvis_ty') || ...
                strcmp(joints{jj},'pelvis_tz')
            %Get the 'force'
            for tt = 0:idData.getDependentColumn([joints{jj},'_force']).size()-1
                ID.(joints{jj})(tt+1,1) = idData.getDependentColumn([joints{jj},'_force']).get(tt);
            end
            clear tt
        else
            %Get the 'moment'
            for tt = 0:idData.getDependentColumn([joints{jj},'_moment']).size()-1
                ID.(joints{jj})(tt+1,1) = idData.getDependentColumn([joints{jj},'_moment']).get(tt);
            end
            clear tt
        end
    end
    clear jj

    %%%%% TODO: consider filter for ID data?????
    
    %Set 'all' structure
    ID.all(:,1) = ID.time;
    for jj = 1:length(joints)
        ID.all(:,jj+1) = ID.(joints{jj});
    end
    clear jj

end