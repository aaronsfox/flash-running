function kinematicsToStates(kinematicsFileName, osimModelFileName,...
    outputFileName, inDegrees, outDegrees)
    
    % Convenience function for converting IK results to a states storage.
    %
    % Input:    kinematicsFileName - file containing kinematic data. Header should only be coordinates name, rather than path to state
    %           osimModelFileName - opensim model filename that corresponds to kinematic data
    %           outputFileName - optional filename to output to (defaults to coordinates.sto)
    %           inDegrees - set to true if kinematics file is in degrees (defaults to True)
    %           outDegrees - set to true if desired output is in degrees (defaults to False)
    
    %Import opensim libraries
    import org.opensim.modeling.*
    
    %Check inputs
    if nargin < 1
        error('Filename for kinematics is required');
    end
    if nargin < 2
        error('OpenSim model filename is required');
    end
    if nargin < 3
        outputFileName = 'coordinates.sto';
    end
    if nargin < 4
        inDegrees = true;
    end
    if nargin < 5
        outDegrees = false;
    end

    %Load in the kinematic data
    kinematicsStorage = Storage(kinematicsFileName);
    
    %Create a copy of the kinematics data to alter the column labels in
    statesStorage = Storage(kinematicsFileName);
    
    %Resample the data points linearly to avoid any later issues with matching
    %time points. Use a time stamp for 250 Hz
    kinematicsStorage.resampleLinear(1/250);
    statesStorage.resampleLinear(1/250);
    
    %Get the column headers for the storage file
    angleNames = kinematicsStorage.getColumnLabels();
    
    %Get the corresponding full paths from the model to rename the
    %angles in the kinematics file
    kinematicModel = Model(osimModelFileName);
    for ii = 0:angleNames.getSize()-1
        currAngle = angleNames.get(ii);
        if ~strcmp(char(currAngle),'time')
            %Get full path to coordinate
            fullPath = [char(kinematicModel.updCoordinateSet().get(currAngle).getAbsolutePathString()),'/value'];
            %Set angle name appropriately using full path
            angleNames.set(ii,fullPath);
        end
    end
    
    %Set the states storage object to have the updated column labels
    statesStorage.setColumnLabels(angleNames);
    
    %Appropriately set output in degrees or radians
    if inDegrees && ~outDegrees
        %Convert degrees values to radians for consistency with the current
        %file label (defaults back to inDegrees=no). Radians seem to work
        %better with the Moco process as well.
        kinematicModel.initSystem();
        kinematicModel.getSimbodyEngine().convertDegreesToRadians(statesStorage);
    elseif inDegrees && outDegrees
        %Change the storage label back to specifying indegrees=yes
        statesStorage.setInDegrees(true);
    elseif ~inDegrees && outDegrees;
        %Convert radians to degrees
        kinematicModel.initSystem();
        kinematicModel.getSimbodyEngine().convertRadiansToDegrees(statesStorage);
        %Reset labeling for degrees
        statesStorage.setInDegrees(true);
    end
    
    %Write the states storage object to file
    statesStorage.print(outputFileName);
    
end