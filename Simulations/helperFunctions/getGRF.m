% This function returns the GRFs given the path to the GRF file
%
% Author: Antoine Falisse
% Date: 12/19/2018
%  
% Modified: Aaron Fox
% Date: 09/11/2020
%
function GRF = getGRF(pathGRF,fileGRF)

    import org.opensim.modeling.*
    
    %Extract ground reaction forces

    %Load grf .mot file
    grfData = TimeSeriesTable([pathGRF,'\',fileGRF]);

    %Get time
    for tt = 0:grfData.getIndependentColumn().size()-1
        GRF.time(tt+1,1) = grfData.getIndependentColumn().get(tt);
    end
    clear tt

    %Set identifiers
    forceIdentifier = {'ground_force_r_v','ground_force_l_v'};
    pointIdentifier = {'ground_force_r_p','ground_force_l_p'};
    momentIdentifier = {'ground_torque_r_','ground_torque_l_'};

    %Set axes and leg variables
    axis = {'x','y','z'};
    leg = {'r','l'};

    %Loop through legs
    for kk = 1:length(leg)
        %Loop through axis
        for aa = 1:length(axis)
            %Get current force data
            currForceLabel = [forceIdentifier{kk},axis{aa}];
            for tt = 0:grfData.getDependentColumn(currForceLabel).size()-1
                GRF.val.(leg{kk})(tt+1,aa) = grfData.getDependentColumn(currForceLabel).get(tt);
            end
            clear tt
            %Get current position data
            currPosLabel = [pointIdentifier{kk},axis{aa}];
            for tt = 0:grfData.getDependentColumn(currPosLabel).size()-1
                GRF.pos.(leg{kk})(tt+1,aa) = grfData.getDependentColumn(currPosLabel).get(tt);
            end
            clear tt
            %Get current moment data
            currMomLabel = [momentIdentifier{kk},axis{aa}];
            for tt = 0:grfData.getDependentColumn(currMomLabel).size()-1
                GRF.Mcop.(leg{kk})(tt+1,aa) = grfData.getDependentColumn(currMomLabel).get(tt);
            end
            clear tt        
        end
        clear aa
    end
    clear kk
    
    %Filter GRF data with 12Hz filter to smooth
    %Filter settings
    order = 4;
    %%%%% TODO: 6Hz vs. 12Hz filter
    cutoff_low = 6;% 12;
    fs = 1/mean(diff(GRF.time(:,1)));
    [af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
    %Loop through legs
    for kk = 1:length(leg)
        %Loop through axis
        for aa = 1:length(axis)
            GRF.val.(leg{kk})(:,aa) = filtfilt(af,bf,GRF.val.(leg{kk})(:,aa));
            GRF.pos.(leg{kk})(:,aa) = filtfilt(af,bf,GRF.pos.(leg{kk})(:,aa));
            GRF.Mcop.(leg{kk})(:,aa) = filtfilt(af,bf,GRF.Mcop.(leg{kk})(:,aa));
        end
        clear aa
    end
    clear kk
    
    %Set data in the 'all' structures
    GRF.val.all = [GRF.time,GRF.val.(leg{1}),GRF.val.(leg{2})];
    GRF.pos.all = [GRF.time,GRF.pos.(leg{1}),GRF.pos.(leg{2})];
    GRF.Mcop.all = [GRF.time,GRF.Mcop.(leg{1}),GRF.Mcop.(leg{2})];
    GRF.Mcop.Y = GRF.Mcop.all(:,[3,6]); 

    %Calculate moments with respect to ground frame origin
    GRF.MorGF.all(:,1) = GRF.time;
    for j = 1:length(leg)
        GRF.MorGF.(leg{j})(:,1) =  GRF.pos.(leg{j})(:,2).*GRF.val.(leg{j})(:,3) - GRF.pos.(leg{j})(:,3).*GRF.val.(leg{j})(:,2);
        GRF.MorGF.(leg{j})(:,2) =  GRF.pos.(leg{j})(:,3).*GRF.val.(leg{j})(:,1) - GRF.pos.(leg{j})(:,1).*GRF.val.(leg{j})(:,3) + GRF.Mcop.(leg{j})(:,2);
        GRF.MorGF.(leg{j})(:,3) =  GRF.pos.(leg{j})(:,1).*GRF.val.(leg{j})(:,2) - GRF.pos.(leg{j})(:,2).*GRF.val.(leg{j})(:,1);
    end
    GRF.MorGF.all(:,2:4) = GRF.MorGF.r;
    GRF.MorGF.all(:,5:7) = GRF.MorGF.l;

end

