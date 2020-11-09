% This function returns the Qs (radian) given the path to the inverse
% kinematics file and the joints of interest.
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
% Modified: Aaron Fox
% Date: 09/11/2020
%
function Qs = getIK(pathIK,fileIK,joints)
    
    import org.opensim.modeling.*

    %Extract joint kinematics
    %Load ik .mot file
    ikData = TimeSeriesTable([pathIK,'\',fileIK]);
    %Get time
    for tt = 0:ikData.getIndependentColumn().size()-1
        Qs.time(tt+1,1) = ikData.getIndependentColumn().get(tt);
    end
    clear tt
    %Get kinematic data
    %Convert to radians as part of this process (except for translations)
    for jj = 1:length(joints)
        %Check for translational coordinate
        if strcmp(joints{jj},'pelvis_tx') || ...
                strcmp(joints{jj},'pelvis_ty') || ...
                strcmp(joints{jj},'pelvis_tz')
            %Don't convert to radians
            for tt = 0:ikData.getDependentColumn(joints{jj}).size()-1
                Qs.(joints{jj})(tt+1,1) = ikData.getDependentColumn(joints{jj}).get(tt);
            end
            clear tt
        else
            %Convert to radians as part of extraction
            for tt = 0:ikData.getDependentColumn(joints{jj}).size()-1
                Qs.(joints{jj})(tt+1,1) = deg2rad(ikData.getDependentColumn(joints{jj}).get(tt));
            end
            clear tt
        end
    end
    clear jj
    %Filter kinematic data with 12Hz filter to match ID data
    %Filter settings
    order = 4;
    %%%%% TODO: 6Hz vs. 12Hz filter
    cutoff_low = 6;% 12;
    fs = 1/mean(diff(Qs.time(:,1)));
    [af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
    %Set time variable in filter, unfilter and colheaders structure
    Qs.allfilt(:,1) = Qs.time;
    Qs.all(:,1) = Qs.time;
    Qs.colheaders{1,1} = 'time';
    %Filter kinematics and store in structure
    %Also use this to store in the 'all' unfiltered structure
    %Also set colheaders variable
    for jj = 1:length(joints)
        Qs.allfilt(:,jj+1) = filtfilt(af,bf,Qs.(joints{jj})); %filtered
        Qs.all(:,jj+1) = Qs.(joints{jj}); %unfiltered
        Qs.colheaders{1,jj+1} = joints{jj}; %colheaders
    end
    clear jj

end
