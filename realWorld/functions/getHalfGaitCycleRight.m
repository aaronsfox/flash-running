%Get half gait cycle for right leg
%(i.e. from right foot strike to left foot strike)
function [startTime, endTime] = getHalfGaitCycleRight(grfFile)
    
    % Convenience function for getting the first half gait cycle based on GRF
    % data. Note that this function is only applicable to the way the current
    % data is structured, but could easily be edited (e.g. for getting a full
    % gait cycle, or left half gait cycle etc.)
    %
    % This function is currently only applicable to going from right contact
    % to another left contact. This could easily be edited with an input 
    % variable asking for limb.
    %
    % Input:    grfFile - .mot file containing GRF time history
    
    %Import opensim libraries
    import org.opensim.modeling.*
    
    %Check input
    if nargin < 1
        error('GRF file is required');
    end
    if ~ischar(grfFile)
        error('Character string of grf file is required');
    end
    if ~strcmpi(grfFile(end-3:end),'.mot')
        error('.mot file extension is required')
    end
        
    %Load the GRF data
    grfTable = TimeSeriesTable(grfFile);
    
    %Convert to more easily readable object for easier manipulation
    %Get number of column labels and data rows
    nLabels = grfTable.getNumColumns();
    nRows = grfTable.getNumRows();
    %Pre-allocate numpy array based on data size
    grfArray = zeros(nRows,nLabels);
    %Loop through labels and rows and get data
    for iLabel = 0:nLabels-1
        %Get current column label
        currLabel = grfTable.getColumnLabel(iLabel);
        %Store column index value if one of the vertical GRF data
        if strcmp(currLabel,'ground_force_r_vy')
            rightGrfCol = iLabel+1;
        elseif strcmp(currLabel,'ground_force_l_vy')
            leftGrfCol = iLabel+1;
        end
        for iRow = 0:nRows-1
            grfArray(iRow+1,iLabel+1) = grfTable.getDependentColumn(currLabel).getElt(0,iRow);
        end
    end
    
    %Create a logical for where the right & left foot is in contact with the ground
    %based on 10N threshold
    %Identify columns for right and left vertical GRF
    rightGrfOn = grfArray(:,rightGrfCol) > 10;
    leftGrfOn = grfArray(:,leftGrfCol) > 10;
    
    %Identify the index where right and left GRF starts
    %Identify where change in boolean value is present
    rightGrfOnInd = find(diff(rightGrfOn) == 1);
    leftGrfOnInd = find(diff(leftGrfOn) == 1);
    
    %Find the first right strike and then get the next left.
    rightStartInd = rightGrfOnInd(1) - 1;
    leftStartInd = leftGrfOnInd(find(leftGrfOnInd > rightStartInd,1))-1;
        
    %Identify the times corresponding to these foot strikes
    startTime = grfTable.getIndependentColumn().get(rightStartInd).doubleValue();
    endTime = grfTable.getIndependentColumn().get(leftStartInd).doubleValue();
    
    %Print outputs
    disp(['Start Time: ',num2str(startTime)]);
    disp(['End Time: ',num2str(endTime)]);
    
end