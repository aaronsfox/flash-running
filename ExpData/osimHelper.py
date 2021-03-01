# Opensim helper functions for processing data

import opensim as osim
import numpy as np
import pandas as pd
import os

# %% addVirtualMarkersStatic

def addVirtualMarkersStatic(staticTRC = None, outputTRC = 'static_withVirtualMarkers.trc'):
    
    # Convenience function for adding virtual markers to static trial.
    #
    # Input:    staticTRC - .trc filename for static trial to add markers to
    #           outputTRC - optional input for filename for output .trc file
    #
    # Hip joint centres are placed according to the method outlined by
    # Harrington et al. (2007), J Biomech, 40: 595-602. 
    #
    # Knee joint centres are placed at the mid point of the femoral epicondyle
    # markers.
    #
    # Pseudo ankle joint centres are placed at the mid points of the malleoli
    # markers.
    #
    # A marker at the mid point of the two metatarsal markers.
    #
    # Markers around the foot and ankle are projected to the floor to assist
    # in aligning the scaling parameters with the axis they are relevant for
    # scaling.
    #
    # Markers at the middle of the set of upper torso and pelvis markers used
    # to help scale torso and pelvis length.
    #
    # Note that the application of this is only relevant to the markerset used
    # with this data. Different labelling of markers will require adaptating
    # this function to make it work.
    
    #Check for input
    if staticTRC is None:
        raise ValueError('Input of staticTRC is required')
    
    #Use the Vec3 TimeSeriesTable to read the Vec3 type data file.
    staticTable = osim.TimeSeriesTableVec3(staticTRC)
    #Convert to more easily readable object for easier manipulation
    #Get number of column labels and data rows
    nLabels = staticTable.getNumColumns()
    nRows = staticTable.getNumRows()
    #Pre-allocate numpy array based on data size
    dataArray = np.zeros((nRows,nLabels*3))
    #Loop through labels and rows and get data
    for iLabel in range(0,nLabels):
        #Get current column label
        currLabel = staticTable.getColumnLabel(iLabel)
        for iRow in range(0,nRows):
            dataArray[iRow,iLabel*3] = staticTable.getDependentColumn(currLabel).getElt(0,iRow).get(0)
            dataArray[iRow,iLabel*3+1] = staticTable.getDependentColumn(currLabel).getElt(0,iRow).get(1)
            dataArray[iRow,iLabel*3+2] = staticTable.getDependentColumn(currLabel).getElt(0,iRow).get(2)
    
    #Convert numpy array to pandas dataframe and add column labels
    #Get column labels
    colLabels = list()
    for iLabel in range(0,nLabels):
        colLabels.append(staticTable.getColumnLabel(iLabel)+'_x')
        colLabels.append(staticTable.getColumnLabel(iLabel)+'_y')
        colLabels.append(staticTable.getColumnLabel(iLabel)+'_z')
    #Convert to dataframe
    static_df = pd.DataFrame(data = dataArray,columns=colLabels)
    
    #Get the pelvis marker data
    #In this step we convert back to the traditional Vicon coordinate system
    #and millimetres. It was easier to do this than mess around with the hip 
    #joint centre calculations
    RASIS = np.zeros((nRows,3))
    RASIS[:,0] = static_df['RASI_z']*1000
    RASIS[:,1] = static_df['RASI_x']*1000
    RASIS[:,2] = static_df['RASI_y']*1000
    LASIS = np.zeros((nRows,3))
    LASIS[:,0] = static_df['LASI_z']*1000
    LASIS[:,1] = static_df['LASI_x']*1000
    LASIS[:,2] = static_df['LASI_y']*1000
    RPSIS = np.zeros((nRows,3))
    RPSIS[:,0] = static_df['RPSI_z']*1000
    RPSIS[:,1] = static_df['RPSI_x']*1000
    RPSIS[:,2] = static_df['RPSI_y']*1000
    LPSIS = np.zeros((nRows,3))
    LPSIS[:,0] = static_df['LPSI_z']*1000
    LPSIS[:,1] = static_df['LPSI_x']*1000
    LPSIS[:,2] = static_df['LPSI_y']*1000
    
    #Calculate hip joint centre at each time step
    #Pre-allocate size
    RHJC = np.zeros((nRows,3))
    LHJC = np.zeros((nRows,3))
    #Loop through sample points
    for t in range(0,nRows):
        #Right handed pelvis reference system definition
        SACRUM = (RPSIS[t,:] + LPSIS[t,:]) / 2
        
        #Global Pelvis Center position
        OP = (LASIS[t,:] + RASIS[t,:]) / 2  
        PROVV = (RASIS[t,:] - SACRUM) / np.linalg.norm(RASIS[t,:] - SACRUM,2)
        IB = (RASIS[t,:] - LASIS[t,:]) / np.linalg.norm(RASIS[t,:] - LASIS[t,:],2)
        KB = np.cross(IB,PROVV)                               
        KB = KB / np.linalg.norm(KB,2)
        JB = np.cross(KB,IB)
        JB = JB / np.linalg.norm(JB,2)
        OB = OP
        
        #rotation+ traslation in homogeneous coordinates (4x4)
        addPelvis = np.array([0,0,0,1])
        pelvis = np.hstack((IB.reshape(3,1),
                            JB.reshape(3,1),
                            KB.reshape(3,1),
                            OB.reshape(3,1)))
        pelvis = np.vstack((pelvis,addPelvis))
        
        #Trasformation into pelvis coordinate system (CS)
        OPB = np.linalg.inv(pelvis) @ np.vstack((OB.reshape(3,1),np.array([1])))    
        PW = np.linalg.norm(RASIS[t,:] - LASIS[t,:])
        PD = np.linalg.norm(SACRUM - OP)
        
        #Harrington formulae (starting from pelvis center)
        diff_ap = -0.24 * PD - 9.9
        diff_v = -0.30 * PW - 10.9
        diff_ml = 0.33 * PW + 7.3    
        
        #vector that must be subtract to OP to obtain hjc in pelvis CS
        vett_diff_pelvis_sx = np.array([-diff_ml,diff_ap,diff_v,1])
        vett_diff_pelvis_dx = np.array([diff_ml,diff_ap,diff_v,1])
        
        #hjc in pelvis CS (4x4)
        rhjc_pelvis = OPB[:,0] + vett_diff_pelvis_dx  
        lhjc_pelvis = OPB[:,0] + vett_diff_pelvis_sx 
        
        #Transformation Local to Global
        RHJC[t,:] = pelvis[0:3,0:3] @ rhjc_pelvis[0:3] + OB
        LHJC[t,:] = pelvis[0:3,0:3] @ lhjc_pelvis[0:3] + OB

    #Add to data frame
    #Convert back to appropriate units and coordinate system here too
    static_df['RHJC_x'] = RHJC[:,1] / 1000
    static_df['RHJC_y'] = RHJC[:,2] / 1000
    static_df['RHJC_z'] = RHJC[:,0] / 1000
    static_df['LHJC_x'] = LHJC[:,1] / 1000
    static_df['LHJC_y'] = LHJC[:,2] / 1000
    static_df['LHJC_z'] = LHJC[:,0] / 1000
    
    #Create mid point markers
    
    #Mid epicondyles
    static_df['RKJC_x'] = (static_df['RLFC_x'] + static_df['RMFC_x']) / 2
    static_df['RKJC_y'] = (static_df['RLFC_y'] + static_df['RMFC_y']) / 2
    static_df['RKJC_z'] = (static_df['RLFC_z'] + static_df['RMFC_z']) / 2
    static_df['LKJC_x'] = (static_df['LLFC_x'] + static_df['LMFC_x']) / 2
    static_df['LKJC_y'] = (static_df['LLFC_y'] + static_df['LMFC_y']) / 2
    static_df['LKJC_z'] = (static_df['LLFC_z'] + static_df['LMFC_z']) / 2
    
    #Mid malleoli
    static_df['RAJC_x'] = (static_df['RLMAL_x'] + static_df['RMMAL_x']) / 2
    static_df['RAJC_y'] = (static_df['RLMAL_y'] + static_df['RMMAL_y']) / 2
    static_df['RAJC_z'] = (static_df['RLMAL_z'] + static_df['RMMAL_z']) / 2
    static_df['LAJC_x'] = (static_df['LLMAL_x'] + static_df['LMMAL_x']) / 2
    static_df['LAJC_y'] = (static_df['LLMAL_y'] + static_df['LMMAL_y']) / 2
    static_df['LAJC_z'] = (static_df['LLMAL_z'] + static_df['LMMAL_z']) / 2
    
    #Mid metatarsal
    static_df['RMT3_x'] = (static_df['RMT1_x'] + static_df['RMT5_x']) / 2
    static_df['RMT3_y'] = (static_df['RMT1_y'] + static_df['RMT5_y']) / 2
    static_df['RMT3_z'] = (static_df['RMT1_z'] + static_df['RMT5_z']) / 2
    static_df['LMT3_x'] = (static_df['LMT1_x'] + static_df['LMT5_x']) / 2
    static_df['LMT3_y'] = (static_df['LMT1_y'] + static_df['LMT5_y']) / 2
    static_df['LMT3_z'] = (static_df['LMT1_z'] + static_df['LMT5_z']) / 2
    
    #Create projections of foot and ankle markers to floor (i.e. y = 0)
    
    #Heel markers
    static_df['fRCAL_x'] = static_df['RCAL_x']
    static_df['fRCAL_y'] = np.zeros((len(static_df),1))
    static_df['fRCAL_z'] = static_df['RCAL_z']
    static_df['fLCAL_x'] = static_df['LCAL_x']
    static_df['fLCAL_y'] = np.zeros((len(static_df),1))
    static_df['fLCAL_z'] = static_df['LCAL_z']
    
    #Ankle joint centres
    static_df['fRAJC_x'] = static_df['RAJC_x']
    static_df['fRAJC_y'] = np.zeros((len(static_df),1))
    static_df['fRAJC_z'] = static_df['RAJC_z']
    static_df['fLAJC_x'] = static_df['LAJC_x']
    static_df['fLAJC_y'] = np.zeros((len(static_df),1))
    static_df['fLAJC_z'] = static_df['LAJC_z']
    
    #Toe markers
    static_df['fRMT1_x'] = static_df['RMT1_x']
    static_df['fRMT1_y'] = np.zeros((len(static_df),1))
    static_df['fRMT1_z'] = static_df['RMT1_z']
    static_df['fLMT1_x'] = static_df['LMT1_x']
    static_df['fLMT1_y'] = np.zeros((len(static_df),1))
    static_df['fLMT1_z'] = static_df['LMT1_z']
    static_df['fRMT5_x'] = static_df['RMT5_x']
    static_df['fRMT5_y'] = np.zeros((len(static_df),1))
    static_df['fRMT5_z'] = static_df['RMT5_z']
    static_df['fLMT5_x'] = static_df['LMT5_x']
    static_df['fLMT5_y'] = np.zeros((len(static_df),1))
    static_df['fLMT5_z'] = static_df['LMT5_z']
    static_df['fRMT3_x'] = static_df['RMT3_x']
    static_df['fRMT3_y'] = np.zeros((len(static_df),1))
    static_df['fRMT3_z'] = static_df['RMT3_z']
    static_df['fLMT3_x'] = static_df['LMT3_x']
    static_df['fLMT3_y'] = np.zeros((len(static_df),1))
    static_df['fLMT3_z'] = static_df['LMT3_z']
    
    #Extra mid points for torso and pelvis markers
    
    #Mid torso
    static_df['MTOR_x'] = (((static_df['RACR_x'] + static_df['LACR_x']) / 2) + ((static_df['C7_x'] + static_df['CLAV_x']) / 2)) / 2
    static_df['MTOR_y'] = (((static_df['RACR_y'] + static_df['LACR_y']) / 2) + ((static_df['C7_y'] + static_df['CLAV_y']) / 2)) / 2
    static_df['MTOR_z'] = (((static_df['RACR_z'] + static_df['LACR_z']) / 2) + ((static_df['C7_z'] + static_df['CLAV_z']) / 2)) / 2
    
    #Mid pelvis
    static_df['MPEL_x'] = (((static_df['RASI_x'] + static_df['LASI_x']) / 2) + ((static_df['RPSI_x'] + static_df['LPSI_x']) / 2)) / 2
    static_df['MPEL_y'] = (((static_df['RASI_y'] + static_df['LASI_y']) / 2) + ((static_df['RPSI_y'] + static_df['LPSI_y']) / 2)) / 2
    static_df['MPEL_z'] = (((static_df['RASI_z'] + static_df['LASI_z']) / 2) + ((static_df['RPSI_z'] + static_df['LPSI_z']) / 2)) / 2
    
    #Build the new time series Vec3 table
    
    #Get the time array
    time = staticTable.getIndependentColumn()
    
    #Set the label names
    #Get the original labels
    markerLabels = list(staticTable.getColumnLabels())
    #Append new marker labels
    markerLabels.append('RHJC')
    markerLabels.append('LHJC')
    markerLabels.append('RKJC')
    markerLabels.append('LKJC')
    markerLabels.append('RAJC')
    markerLabels.append('LAJC')
    markerLabels.append('RMT3')
    markerLabels.append('LMT3')
    markerLabels.append('fRCAL')
    markerLabels.append('fLCAL')
    markerLabels.append('fRAJC')
    markerLabels.append('fLAJC')
    markerLabels.append('fRMT1')
    markerLabels.append('fLMT1')
    markerLabels.append('fRMT5')
    markerLabels.append('fLMT5')
    markerLabels.append('fRMT3')
    markerLabels.append('fLMT3')
    markerLabels.append('MTOR')
    markerLabels.append('MPEL')
    
    #Initialise time series table of vec3 format
    newTable = osim.TimeSeriesTableVec3()
    
    #Set labels in table
    newLabels = osim.StdVectorString()
    for ii in range(0,len(markerLabels)):
        newLabels.append(markerLabels[ii])
    newTable.setColumnLabels(newLabels)
    
    #Build the time series table
    #Loop through rows and allocate data
    for iRow in range(0,nRows):
        #Create a blank opensim row vector
        row = osim.RowVectorVec3(len(markerLabels))
        #Create and fill row data
        for iCol in range(0,len(markerLabels)):
            #Identify the dataframe labels for the current marker label
            xLabel = markerLabels[iCol]+'_x'
            yLabel = markerLabels[iCol]+'_y'
            zLabel = markerLabels[iCol]+'_z'
            #Set data in row
            row.getElt(0,iCol).set(0,static_df[xLabel][iRow])
            row.getElt(0,iCol).set(1,static_df[yLabel][iRow])
            row.getElt(0,iCol).set(2,static_df[zLabel][iRow])
        #Append row to table
        newTable.appendRow(iRow,row)
        #Set time value
        newTable.setIndependentValueAtIndex(iRow,time[iRow])
            
    #Update the meta data in the new table based on the original
    
    #Get meta data keys
    metaKeys = list(staticTable.getTableMetaDataKeys())
    
    #Loop through keys and append to new table
    for ii in range(0,len(metaKeys)):
        #Get current meta data string
        currMeta = metaKeys[ii]
        #Append with relevant value
        newTable.addTableMetaDataString(currMeta,staticTable.getTableMetaDataAsString(currMeta))
        
    #Write the new table to file
    osim.TRCFileAdapter().write(newTable,outputTRC)
    
    #Print confirm message
    print('Virtual markers added to new .trc file: '+outputTRC)

# %% getHalfGaitCycle

def getHalfGaitCycle(grfFile = None):
    
    # Convenience function for getting the middle half gait cycle based on GRF
    # data. Note that this function is only applicable to the way the current
    # data is structured, but could easily be edited (e.g. for getting a full
    # gait cycle, or left to right foot strike etc.)
    #
    # Input:    grfFile - .mot file containing GRF time history
    
    #Check for input
    if grfFile is None:
        raise ValueError('A GRF file in .mot format is required')
        
    #Load the GRF data
    grfTable = osim.TimeSeriesTable(grfFile)
    
    #Convert to more easily readable object for easier manipulation
    #Get number of column labels and data rows
    nLabels = grfTable.getNumColumns()
    nRows = grfTable.getNumRows()
    #Pre-allocate numpy array based on data size
    grfArray = np.zeros((nRows,nLabels))
    #Loop through labels and rows and get data
    for iLabel in range(0,nLabels):
        #Get current column label
        currLabel = grfTable.getColumnLabel(iLabel)
        #Store column index value if one of the vertical GRF data
        if currLabel == 'ground_force_r_vy':
            rightGrfCol = int(iLabel)
        elif  currLabel == 'ground_force_l_vy':
            leftGrfCol = int(iLabel)
        for iRow in range(0,nRows):
            grfArray[iRow,iLabel] = grfTable.getDependentColumn(currLabel).getElt(0,iRow)
    
    #Create a logical for where the right & left foot is in contact with the ground
    #based on 10N threshold
    #Identify columns for right and left vertical GRF
    rightGrfOn = grfArray[:,rightGrfCol] > 10
    leftGrfOn = grfArray[:,leftGrfCol] > 10
    
    #Identify the index where right and left GRF starts
    #Identify where change in boolean value is present
    rightGrfDiff = np.where(rightGrfOn[:-1] != rightGrfOn[1:])[0]
    leftGrfDiff = np.where(leftGrfOn[:-1] != leftGrfOn[1:])[0]
    #Identify where the index one above the change is true for foot strikes
    rightGrfOnInd = list()
    for gg in range(0,len(rightGrfDiff)):
        if rightGrfOn[rightGrfDiff[gg]+1]:
            rightGrfOnInd.append(rightGrfDiff[gg]+1)
    leftGrfOnInd = list()
    for gg in range(0,len(leftGrfDiff)):
        if leftGrfOn[leftGrfDiff[gg]+1]:
            leftGrfOnInd.append(leftGrfDiff[gg]+1)
    
    #Find the middle right foot strike
    rightStartInd = rightGrfOnInd[int(len(rightGrfOnInd)/2)] - 1
    
    #Find the next left foot strike after the right foot strike
    leftStartInd = leftGrfOnInd[np.where(leftGrfOnInd > rightStartInd)[0][0]] - 1
    
    #Identify the times corresponding to these foot strikes
    startTime = list(grfTable.getIndependentColumn())[rightStartInd]
    endTime = list(grfTable.getIndependentColumn())[leftStartInd]
    
    #Print outputs
    print('Start Time: '+str(startTime))
    print('End Time: '+str(endTime))
    
    return startTime,endTime

# %% getFullGaitCycle

def getFullGaitCycle(grfFile = None):
    
    # Convenience function for getting the first full gait cycle based on GRF
    # data. Note that this function is only applicable to the way the current
    # data is structured, but could easily be edited (e.g. for getting a full
    # gait cycle, or right to right foot strike etc.)
    #
    # This function is currently only applicable to going from right contact
    # to another right contact. This could easily be edited with an input 
    # variable asking for limb.
    #
    # Input:    grfFile - .mot file containing GRF time history
    
    #Check for input
    if grfFile is None:
        raise ValueError('A GRF file in .mot format is required')
        
    #Load the GRF data
    grfTable = osim.TimeSeriesTable(grfFile)
    
    #Convert to more easily readable object for easier manipulation
    #Get number of column labels and data rows
    nLabels = grfTable.getNumColumns()
    nRows = grfTable.getNumRows()
    #Pre-allocate numpy array based on data size
    grfArray = np.zeros((nRows,nLabels))
    #Loop through labels and rows and get data
    for iLabel in range(0,nLabels):
        #Get current column label
        currLabel = grfTable.getColumnLabel(iLabel)
        #Store column index value if one of the vertical GRF data
        if currLabel == 'ground_force_r_vy':
            rightGrfCol = int(iLabel)
        elif  currLabel == 'ground_force_l_vy':
            leftGrfCol = int(iLabel)
        for iRow in range(0,nRows):
            grfArray[iRow,iLabel] = grfTable.getDependentColumn(currLabel).getElt(0,iRow)
    
    #Create a logical for where the right & left foot is in contact with the ground
    #based on 10N threshold
    #Identify columns for right and left vertical GRF
    rightGrfOn = grfArray[:,rightGrfCol] > 10
    leftGrfOn = grfArray[:,leftGrfCol] > 10
    
    #Identify the index where right and left GRF starts
    #Identify where change in boolean value is present
    rightGrfDiff = np.where(rightGrfOn[:-1] != rightGrfOn[1:])[0]
    leftGrfDiff = np.where(leftGrfOn[:-1] != leftGrfOn[1:])[0]
    #Identify where the index one above the change is true for foot strikes
    rightGrfOnInd = list()
    for gg in range(0,len(rightGrfDiff)):
        if rightGrfOn[rightGrfDiff[gg]+1]:
            rightGrfOnInd.append(rightGrfDiff[gg]+1)
    leftGrfOnInd = list()
    for gg in range(0,len(leftGrfDiff)):
        if leftGrfOn[leftGrfDiff[gg]+1]:
            leftGrfOnInd.append(leftGrfDiff[gg]+1)
    
    #Instead of finding the middle foot strike like in the half gait cycle,
    #find the first and then get the next. With sprint trials there are generally
    #only two anyway
    rightStartInd = rightGrfOnInd[0] - 1
    rightStopInd = rightGrfOnInd[1] - 1
    
    #Identify the times corresponding to these foot strikes
    startTime = list(grfTable.getIndependentColumn())[rightStartInd]
    endTime = list(grfTable.getIndependentColumn())[rightStopInd]
    
    #Print outputs
    print('Start Time: '+str(startTime))
    print('End Time: '+str(endTime))
    
    return startTime,endTime

# %% addVirtualMarkersDynamic

def addVirtualMarkersDynamic(staticTRC = None, dynamicTRC = None,
                             outputTRC = 'dynamic_withVirtualMarkers.trc'):
    
    # Convenience function for adding virtual markers to dynamic trial.
    #
    # Input:    staticTRC - .trc filename for static trial with virtual markers added
    #           dynamicTRC - .trc filename for dynamic trial to add markers to
    #           outputTRC - optional input for filename for output .trc file
    #
    # Hip joint centres are placed according to the method outlined by
    # Harrington et al. (2007), J Biomech, 40: 595-602. 
    #
    # A 3d trilateration function from https://github.com/akshayb6/trilateration-in-3d
    # is used to create virtual markers at positions on the scaled model based 
    # on distances and positions of four other markers.
    #
    # Gaps in marker data come across as zeros in the .trc file, so any calculations
    # over where missing markers are present will be wrong. These areas should
    # not be used in the actual trial assessment without appropriate filling.
    
    #Check for input
    if staticTRC is None:
        raise ValueError('Input of staticTRC is required')
    if dynamicTRC is None:
        raise ValueError('Input of dynamicTRC is required')
    
    #Use the Vec3 TimeSeriesTable to read the Vec3 type data file.
    dynamicTable = osim.TimeSeriesTableVec3(dynamicTRC)
    
    #Convert to more easily readable object for easier manipulation
    #Get number of column labels and data rows
    nLabels = dynamicTable.getNumColumns()
    nRows = dynamicTable.getNumRows()
    #Pre-allocate numpy array based on data size
    dataArray = np.zeros((nRows,nLabels*3))
    #Loop through labels and rows and get data
    for iLabel in range(0,nLabels):
        #Get current column label
        currLabel = dynamicTable.getColumnLabel(iLabel)
        for iRow in range(0,nRows):
            dataArray[iRow,iLabel*3] = dynamicTable.getDependentColumn(currLabel).getElt(0,iRow).get(0)
            dataArray[iRow,iLabel*3+1] = dynamicTable.getDependentColumn(currLabel).getElt(0,iRow).get(1)
            dataArray[iRow,iLabel*3+2] = dynamicTable.getDependentColumn(currLabel).getElt(0,iRow).get(2)
    
    #Convert numpy array to pandas dataframe and add column labels
    #Get column labels
    colLabels = list()
    for iLabel in range(0,nLabels):
        colLabels.append(dynamicTable.getColumnLabel(iLabel)+'_x')
        colLabels.append(dynamicTable.getColumnLabel(iLabel)+'_y')
        colLabels.append(dynamicTable.getColumnLabel(iLabel)+'_z')
    #Convert to dataframe
    dynamic_df = pd.DataFrame(data = dataArray,columns=colLabels)
    
    #Get the pelvis marker data
    #In this step we convert back to the traditional Vicon coordinate system
    #and millimetres. It was easier to do this than mess around with the hip 
    #joint centre calculations
    RASIS = np.zeros((nRows,3))
    RASIS[:,0] = dynamic_df['RASI_z']*1000
    RASIS[:,1] = dynamic_df['RASI_x']*1000
    RASIS[:,2] = dynamic_df['RASI_y']*1000
    LASIS = np.zeros((nRows,3))
    LASIS[:,0] = dynamic_df['LASI_z']*1000
    LASIS[:,1] = dynamic_df['LASI_x']*1000
    LASIS[:,2] = dynamic_df['LASI_y']*1000
    RPSIS = np.zeros((nRows,3))
    RPSIS[:,0] = dynamic_df['RPSI_z']*1000
    RPSIS[:,1] = dynamic_df['RPSI_x']*1000
    RPSIS[:,2] = dynamic_df['RPSI_y']*1000
    LPSIS = np.zeros((nRows,3))
    LPSIS[:,0] = dynamic_df['LPSI_z']*1000
    LPSIS[:,1] = dynamic_df['LPSI_x']*1000
    LPSIS[:,2] = dynamic_df['LPSI_y']*1000
    
    #Calculate hip joint centre at each time step
    #Pre-allocate size
    RHJC = np.zeros((nRows,3))
    LHJC = np.zeros((nRows,3))
    #Loop through sample points
    for t in range(0,nRows):
        #Right handed pelvis reference system definition
        SACRUM = (RPSIS[t,:] + LPSIS[t,:]) / 2
        
        #Global Pelvis Center position
        OP = (LASIS[t,:] + RASIS[t,:]) / 2  
        PROVV = (RASIS[t,:] - SACRUM) / np.linalg.norm(RASIS[t,:] - SACRUM,2)
        IB = (RASIS[t,:] - LASIS[t,:]) / np.linalg.norm(RASIS[t,:] - LASIS[t,:],2)
        KB = np.cross(IB,PROVV)                               
        KB = KB / np.linalg.norm(KB,2)
        JB = np.cross(KB,IB)
        JB = JB / np.linalg.norm(JB,2)
        OB = OP
        
        #rotation+ traslation in homogeneous coordinates (4x4)
        addPelvis = np.array([0,0,0,1])
        pelvis = np.hstack((IB.reshape(3,1),
                            JB.reshape(3,1),
                            KB.reshape(3,1),
                            OB.reshape(3,1)))
        pelvis = np.vstack((pelvis,addPelvis))
        
        #Trasformation into pelvis coordinate system (CS)
        OPB = np.linalg.inv(pelvis) @ np.vstack((OB.reshape(3,1),np.array([1])))    
        PW = np.linalg.norm(RASIS[t,:] - LASIS[t,:])
        PD = np.linalg.norm(SACRUM - OP)
        
        #Harrington formulae (starting from pelvis center)
        diff_ap = -0.24 * PD - 9.9
        diff_v = -0.30 * PW - 10.9
        diff_ml = 0.33 * PW + 7.3    
        
        #vector that must be subtract to OP to obtain hjc in pelvis CS
        vett_diff_pelvis_sx = np.array([-diff_ml,diff_ap,diff_v,1])
        vett_diff_pelvis_dx = np.array([diff_ml,diff_ap,diff_v,1])
        
        #hjc in pelvis CS (4x4)
        rhjc_pelvis = OPB[:,0] + vett_diff_pelvis_dx  
        lhjc_pelvis = OPB[:,0] + vett_diff_pelvis_sx 
        
        #Transformation Local to Global
        RHJC[t,:] = pelvis[0:3,0:3] @ rhjc_pelvis[0:3] + OB
        LHJC[t,:] = pelvis[0:3,0:3] @ lhjc_pelvis[0:3] + OB

    #Add to data frame
    #Convert back to appropriate units and coordinate system here too
    dynamic_df['RHJC_x'] = RHJC[:,1] / 1000
    dynamic_df['RHJC_y'] = RHJC[:,2] / 1000
    dynamic_df['RHJC_z'] = RHJC[:,0] / 1000
    dynamic_df['LHJC_x'] = LHJC[:,1] / 1000
    dynamic_df['LHJC_y'] = LHJC[:,2] / 1000
    dynamic_df['LHJC_z'] = LHJC[:,0] / 1000
    
    #Create mid point markers
    
    #Mid torso
    dynamic_df['MTOR_x'] = (((dynamic_df['RACR_x'] + dynamic_df['LACR_x']) / 2) + ((dynamic_df['C7_x'] + dynamic_df['CLAV_x']) / 2)) / 2
    dynamic_df['MTOR_y'] = (((dynamic_df['RACR_y'] + dynamic_df['LACR_y']) / 2) + ((dynamic_df['C7_y'] + dynamic_df['CLAV_y']) / 2)) / 2
    dynamic_df['MTOR_z'] = (((dynamic_df['RACR_z'] + dynamic_df['LACR_z']) / 2) + ((dynamic_df['C7_z'] + dynamic_df['CLAV_z']) / 2)) / 2
    
    #Mid pelvis
    dynamic_df['MPEL_x'] = (((dynamic_df['RASI_x'] + dynamic_df['LASI_x']) / 2) + ((dynamic_df['RPSI_x'] + dynamic_df['LPSI_x']) / 2)) / 2
    dynamic_df['MPEL_y'] = (((dynamic_df['RASI_y'] + dynamic_df['LASI_y']) / 2) + ((dynamic_df['RPSI_y'] + dynamic_df['LPSI_y']) / 2)) / 2
    dynamic_df['MPEL_z'] = (((dynamic_df['RASI_z'] + dynamic_df['LASI_z']) / 2) + ((dynamic_df['RPSI_z'] + dynamic_df['LPSI_z']) / 2)) / 2
    
    #Create virtual markers using 3d trilateration
    
    #Define the 3d trilateration function
    #Trilateration calculations (from: https://github.com/akshayb6/trilateration-in-3d)
    def trilateration3D(p1 = None, p2 = None, p3 = None, p4 = None,
                        d1 = None, d2 = None, d3 = None, d4 = None):
        
        # Inputs
        # p1, p2, p3 and p4 - 3 column x nrow array of the 3d marker points
        # d1, d2, d3 and d4 - distances to the virtual marker in the static trial from the corresponding points
        
        #Check inputs
        if p1 is None or p2 is None or p3 is None or p4 is None:
            raise ValueError('All 4 points need to be specified as nd array')
        if d1 is None or d2 is None or d3 is None or d4 is None:
            raise ValueError('All 4 distances need to be specified as floats')
        
        #Allocate numpy array for results
        virtualMarker = np.zeros((len(p1),3))
        
        #Run calculations
        #Loop through numpy arrays
        for pp in range(0,len(p1)):
            e_x = (p2[pp]-p1[pp])/np.linalg.norm(p2[pp]-p1[pp])
            i = np.dot(e_x,(p3[pp]-p1[pp]))
            e_y = (p3[pp]-p1[pp]-(i*e_x))/(np.linalg.norm(p3[pp]-p1[pp]-(i*e_x)))
            e_z = np.cross(e_x,e_y)
            d = np.linalg.norm(p2[pp]-p1[pp])
            j = np.dot(e_y,(p3[pp]-p1[pp]))
            x = ((d1**2)-(d2**2)+(d**2))/(2*d)
            y = (((d1**2)-(d3**2)+(i**2)+(j**2))/(2*j))-((i/j)*(x))
            z1 = np.sqrt(d1**2-x**2-y**2)
            z2 = np.sqrt(d1**2-x**2-y**2)*(-1)
            ans1 = p1[pp]+(x*e_x)+(y*e_y)+(z1*e_z)
            ans2 = p1[pp]+(x*e_x)+(y*e_y)+(z2*e_z)
            dist1 = np.linalg.norm(p4[pp]-ans1)
            dist2 = np.linalg.norm(p4[pp]-ans2)
            if np.abs(d4-dist1) < np.abs(d4-dist2):
                virtualMarker[pp,:] = ans1
            else: 
                virtualMarker[pp,:] = ans2
        #Return calculated virtual marker from function
        return virtualMarker            
    
    #Load in the static marker file to identify average marker distances
    #Use the Vec3 TimeSeriesTable to read the Vec3 type data file.
    staticTable = osim.TimeSeriesTableVec3(staticTRC)
    #Convert to more easily readable object for easier manipulation
    #Get number of column labels and data rows
    nLabelsStatic = staticTable.getNumColumns()
    nRowsStatic = staticTable.getNumRows()
    #Pre-allocate numpy array based on data size
    dataArrayStatic = np.zeros((nRowsStatic,nLabelsStatic*3))
    #Loop through labels and rows and get data
    for iLabel in range(0,nLabelsStatic):
        #Get current column label
        currLabel = staticTable.getColumnLabel(iLabel)
        for iRow in range(0,nRowsStatic):
            dataArrayStatic[iRow,iLabel*3] = staticTable.getDependentColumn(currLabel).getElt(0,iRow).get(0)
            dataArrayStatic[iRow,iLabel*3+1] = staticTable.getDependentColumn(currLabel).getElt(0,iRow).get(1)
            dataArrayStatic[iRow,iLabel*3+2] = staticTable.getDependentColumn(currLabel).getElt(0,iRow).get(2)
    
    #Convert numpy array to pandas dataframe and add column labels
    #Get column labels
    colLabelsStatic = list()
    for iLabel in range(0,nLabelsStatic):
        colLabelsStatic.append(staticTable.getColumnLabel(iLabel)+'_x')
        colLabelsStatic.append(staticTable.getColumnLabel(iLabel)+'_y')
        colLabelsStatic.append(staticTable.getColumnLabel(iLabel)+'_z')
    #Convert to dataframe
    static_df = pd.DataFrame(data = dataArrayStatic,columns=colLabelsStatic)
    
    #Right and left knee joint centres
    #This will use the location and distances of the hip joint centre, and thigh
    #cluster markers
    
    #Calculate the mean distances for the four markers to the added knee joint centre
    
    # #Right limb
    # #Calculate distances from relevant markers
    # RHJC_RKJC_d = np.sqrt(np.square((np.mean(static_df['RKJC_x']) - np.mean(static_df['RHJC_x']))) + 
    #                       np.square((np.mean(static_df['RKJC_y']) - np.mean(static_df['RHJC_y']))) + 
    #                       np.square((np.mean(static_df['RKJC_z']) - np.mean(static_df['RHJC_z']))))
    # RTH1_RKJC_d = np.sqrt(np.square((np.mean(static_df['RKJC_x']) - np.mean(static_df['RTH1_x']))) + 
    #                       np.square((np.mean(static_df['RKJC_y']) - np.mean(static_df['RTH1_y']))) + 
    #                       np.square((np.mean(static_df['RKJC_z']) - np.mean(static_df['RTH1_z']))))
    # RTH2_RKJC_d = np.sqrt(np.square((np.mean(static_df['RKJC_x']) - np.mean(static_df['RTH2_x']))) + 
    #                       np.square((np.mean(static_df['RKJC_y']) - np.mean(static_df['RTH2_y']))) + 
    #                       np.square((np.mean(static_df['RKJC_z']) - np.mean(static_df['RTH2_z']))))
    # RTH3_RKJC_d = np.sqrt(np.square((np.mean(static_df['RKJC_x']) - np.mean(static_df['RTH3_x']))) + 
    #                       np.square((np.mean(static_df['RKJC_y']) - np.mean(static_df['RTH3_y']))) + 
    #                       np.square((np.mean(static_df['RKJC_z']) - np.mean(static_df['RTH3_z']))))
    # #Run trilateration function
    # RKJC = trilateration3D(p1 = dynamic_df[['RHJC_x','RHJC_y','RHJC_z']].to_numpy(),
    #                        p2 = dynamic_df[['RTH1_x','RTH1_y','RTH1_z']].to_numpy(),
    #                        p3 = dynamic_df[['RTH2_x','RTH2_y','RTH2_z']].to_numpy(),
    #                        p4 = dynamic_df[['RTH3_x','RTH3_y','RTH3_z']].to_numpy(),
    #                        d1 = RHJC_RKJC_d, d2 = RTH1_RKJC_d, d3 = RTH2_RKJC_d, d4 = RTH3_RKJC_d)
    # #Allocate to dataframe
    # dynamic_df['RKJC_x'] = RKJC[:,0]
    # dynamic_df['RKJC_y'] = RKJC[:,1]
    # dynamic_df['RKJC_z'] = RKJC[:,2]
    
    # #Left limb
    # #Calculate distances from relevant markers
    # LHJC_LKJC_d = np.sqrt(np.square((np.mean(static_df['LKJC_x']) - np.mean(static_df['LHJC_x']))) + 
    #                       np.square((np.mean(static_df['LKJC_y']) - np.mean(static_df['LHJC_y']))) + 
    #                       np.square((np.mean(static_df['LKJC_z']) - np.mean(static_df['LHJC_z']))))
    # LTH1_LKJC_d = np.sqrt(np.square((np.mean(static_df['LKJC_x']) - np.mean(static_df['LTH1_x']))) + 
    #                       np.square((np.mean(static_df['LKJC_y']) - np.mean(static_df['LTH1_y']))) + 
    #                       np.square((np.mean(static_df['LKJC_z']) - np.mean(static_df['LTH1_z']))))
    # LTH2_LKJC_d = np.sqrt(np.square((np.mean(static_df['LKJC_x']) - np.mean(static_df['LTH2_x']))) + 
    #                       np.square((np.mean(static_df['LKJC_y']) - np.mean(static_df['LTH2_y']))) + 
    #                       np.square((np.mean(static_df['LKJC_z']) - np.mean(static_df['LTH2_z']))))
    # LTH3_LKJC_d = np.sqrt(np.square((np.mean(static_df['LKJC_x']) - np.mean(static_df['LTH3_x']))) + 
    #                       np.square((np.mean(static_df['LKJC_y']) - np.mean(static_df['LTH3_y']))) + 
    #                       np.square((np.mean(static_df['LKJC_z']) - np.mean(static_df['LTH3_z']))))
    # #Run trilateration function
    # LKJC = trilateration3D(p1 = dynamic_df[['LHJC_x','LHJC_y','LHJC_z']].to_numpy(),
    #                        p2 = dynamic_df[['LTH1_x','LTH1_y','LTH1_z']].to_numpy(),
    #                        p3 = dynamic_df[['LTH2_x','LTH2_y','LTH2_z']].to_numpy(),
    #                        p4 = dynamic_df[['LTH3_x','LTH3_y','LTH3_z']].to_numpy(),
    #                        d1 = LHJC_LKJC_d, d2 = LTH1_LKJC_d, d3 = LTH2_LKJC_d, d4 = LTH3_LKJC_d)
    # #Allocate to dataframe
    # dynamic_df['LKJC_x'] = LKJC[:,0]
    # dynamic_df['LKJC_y'] = LKJC[:,1]
    # dynamic_df['LKJC_z'] = LKJC[:,2]
    
    ##### NOTE: trialteration doesn't work that great, markers move around...

    #Build the new time series Vec3 table
    
    #Get the time array
    time = dynamicTable.getIndependentColumn()
    
    #Set the label names
    #Get the original labels
    markerLabels = list(dynamicTable.getColumnLabels())
    #Append new marker labels
    markerLabels.append('RHJC')
    markerLabels.append('LHJC')
    # markerLabels.append('RKJC')
    # markerLabels.append('LKJC')
    markerLabels.append('MTOR')
    markerLabels.append('MPEL')
    
    #Initialise time series table of vec3 format
    newTable = osim.TimeSeriesTableVec3()
    
    #Set labels in table
    newLabels = osim.StdVectorString()
    for ii in range(0,len(markerLabels)):
        newLabels.append(markerLabels[ii])
    newTable.setColumnLabels(newLabels)
    
    #Build the time series table
    #Loop through rows and allocate data
    for iRow in range(0,nRows):
        #Create a blank opensim row vector
        row = osim.RowVectorVec3(len(markerLabels))
        #Create and fill row data
        for iCol in range(0,len(markerLabels)):
            #Identify the dataframe labels for the current marker label
            xLabel = markerLabels[iCol]+'_x'
            yLabel = markerLabels[iCol]+'_y'
            zLabel = markerLabels[iCol]+'_z'
            #Set data in row
            row.getElt(0,iCol).set(0,dynamic_df[xLabel][iRow])
            row.getElt(0,iCol).set(1,dynamic_df[yLabel][iRow])
            row.getElt(0,iCol).set(2,dynamic_df[zLabel][iRow])
        #Append row to table
        newTable.appendRow(iRow,row)
        #Set time value
        newTable.setIndependentValueAtIndex(iRow,time[iRow])
            
    #Update the meta data in the new table based on the original
    
    #Get meta data keys
    metaKeys = list(dynamicTable.getTableMetaDataKeys())
    
    #Loop through keys and append to new table
    for ii in range(0,len(metaKeys)):
        #Get current meta data string
        currMeta = metaKeys[ii]
        #Append with relevant value
        newTable.addTableMetaDataString(currMeta,dynamicTable.getTableMetaDataAsString(currMeta))
        
    #Write the new table to file
    osim.TRCFileAdapter().write(newTable,outputTRC)
    
    #Print confirm message
    print('Virtual markers added to new .trc file: '+outputTRC)
       
# %% scaleOptimalForceSubjectSpecific
  
def scaleOptimalForceSubjectSpecific(genericModelFileName = None, scaledModelFileName = None,
                                     genericHeight = 1.700, scaledHeight = None,
                                     outputModelFileName = 'strengthScaled.osim'):
    
    # Convenience function for scaling muscle strength based on height/mass relationship.
    # The scaling is based on the Handsfield et al. equations, and this function
    # is adapted from that which comes with the Rajagopal et al. model
    #
    # Input:    genericModelFileName - .osim filename of the original generic model used
    #           scaledModelFileName - .osim filename of the scaled model
    #           genericHeight - height (in m) of the generic model (defaults to Rajagopal model height)
    #           scaledHeight - height (in m) of the participant used to scale model
    
    #Check for input
    if genericModelFileName is None:
        raise ValueError('Filename for generic model is required')
    if scaledModelFileName is None:
        raise ValueError('Filename for scaled model is required')
    if scaledHeight is None:
        raise ValueError('Height for participant (in m) is required')
    
    #Load in models
    genericModel = osim.Model(genericModelFileName)
    scaledModel = osim.Model(scaledModelFileName)
    
    #Get total mass of the two models
    #Generic model
    genericMass = 0
    allBodies = genericModel.getBodySet()
    for ii in range(0,allBodies.getSize()):
        genericMass = genericMass + allBodies.get(ii).getMass()
    #Scaled model
    scaledMass = 0
    allBodies = scaledModel.getBodySet()
    for ii in range(0,allBodies.getSize()):
        scaledMass = scaledMass + allBodies.get(ii).getMass()
    
    #Get muscle volume totals based on mass and heights
    genericVtotal = 47.05 * genericMass * genericHeight + 1289.6
    scaledVtotal = 47.05 * scaledMass * scaledHeight + 1289.6
    
    #Get the muscles
    genericAllMuscles = genericModel.getMuscles()
    scaledAllMuscles = scaledModel.getMuscles()
    
    #Loop through all muscles and scale according to volume and muscle parameters
    for ii in range(0,genericAllMuscles.getSize()):
        #Get current muscle name
        #Using name avoids any issues with differing muscle order, but all names
        #must be the same
        currMuscle = genericAllMuscles.get(ii).getName()
        #Get muscle from each model
        genericCurrentMuscle = genericAllMuscles.get(currMuscle)
        scaledCurrentMuscle = scaledAllMuscles.get(currMuscle)
        #Get optimal fibre length for each muscle
        genericL0 = genericCurrentMuscle.getOptimalFiberLength()
        scaledL0 = scaledCurrentMuscle.getOptimalFiberLength()
        #Set force scale factor
        forceScaleFactor = (scaledVtotal / genericVtotal) / (scaledL0 / genericL0)
        #Scale current muscle strength
        scaledCurrentMuscle.setMaxIsometricForce(forceScaleFactor * scaledCurrentMuscle.getMaxIsometricForce())
        
    #Print out scaled muscle model
    scaledModel.finalizeConnections()
    scaledModel.printToXML(outputModelFileName)
    
# %% kinematicsToStates

def kinematicsToStates(kinematicsFileName = None, osimModelFileName = None,
                       outputFileName = 'coordinates.sto',
                       inDegrees = True, outDegrees = False):
    
    # Convenience function for converting IK results to a states storage.
    #
    # Input:    kinematicsFileName - file containing kinematic data. Header should only be coordinates name, rather than path to state
    #           osimModelFileName - opensim model filename that corresponds to kinematic data
    #           outputFileName - optional filename to output to (defaults to coordinates.sto)
    #           inDegrees - set to true if kinematics file is in degrees (defaults to True)
    #           outDegrees - set to true if desired output is in degrees (defaults to False)

    if kinematicsFileName is None:
        raise ValueError('Filename for kinematics is required')
    if osimModelFileName is None:
        raise ValueError('OpenSim model filename is required')

    #Load in the kinematic data
    kinematicsStorage = osim.Storage(kinematicsFileName)
    
    #Create a copy of the kinematics data to alter the column labels in
    statesStorage = osim.Storage(kinematicsFileName)
    
    #Resample the data points linearly to avoid any later issues with matching
    #time points. Use a time stamp for 250 Hz
    kinematicsStorage.resampleLinear(1/250)
    statesStorage.resampleLinear(1/250)
    
    #Get the column headers for the storage file
    angleNames = kinematicsStorage.getColumnLabels()
    
    #Get the corresponding full paths from the model to rename the
    #angles in the kinematics file
    kinematicModel = osim.Model(osimModelFileName)
    for ii in range(0,angleNames.getSize()):
        currAngle = angleNames.get(ii)
        if currAngle != 'time':
            #Get full path to coordinate
            fullPath = kinematicModel.updCoordinateSet().get(currAngle).getAbsolutePathString()+'/value'
            #Set angle name appropriately using full path
            angleNames.set(ii,fullPath)
    
    #Set the states storage object to have the updated column labels
    statesStorage.setColumnLabels(angleNames)
    
    #Appropriately set output in degrees or radians
    if inDegrees and not outDegrees:
        #Convert degrees values to radians for consistency with the current
        #file label (defaults back to inDegrees=no). Radians seem to work
        #better with the Moco process as well.
        kinematicModel.initSystem()
        kinematicModel.getSimbodyEngine().convertDegreesToRadians(statesStorage)
    elif inDegrees and outDegrees:
        #Change the storage label back to specifying indegrees=yes
        statesStorage.setInDegrees(True)
    elif not inDegrees and outDegrees:
        #Convert radians to degrees
        kinematicModel.initSystem()
        kinematicModel.getSimbodyEngine().convertRadiansToDegrees(statesStorage)
        #Reset labeling for degrees
        statesStorage.setInDegrees(True)
    
    #Write the states storage object to file
    statesStorage.printToXML(outputFileName)

# %% statesTo2D

def statesTo2D(statesFileName = None, outputFileName = 'coordinates.sto'):
    
    # Convenience function for converting the created states from IK to match
    # up to the 2D model. Note that this function requires the input states file
    # to be in radians, as there is no process for converting these.
    #
    # Input:    statesFileName - file containing states data from kinematic conversion
    #           outputFileName - optional filename to output to (defaults to coordinates.sto)

    if statesFileName is None:
        raise ValueError('Filename for states is required')

    #Load in the kinematic data
    origStorage = osim.TimeSeriesTable(statesFileName)
    
    #Set column labels to remove
    removeLabels = ['/jointset/ground_pelvis/pelvis_list/value',
                    '/jointset/ground_pelvis/pelvis_rotation/value',
                    '/jointset/ground_pelvis/pelvis_tz/value',
                    '/jointset/hip_r/hip_adduction_r/value',
                    '/jointset/hip_r/hip_rotation_r/value',
                    '/jointset/patellofemoral_r/knee_angle_r_beta/value',
                    '/jointset/hip_l/hip_adduction_l/value',
                    '/jointset/hip_l/hip_rotation_l/value',
                    '/jointset/patellofemoral_l/knee_angle_l_beta/value',
                    '/jointset/back/lumbar_bending/value',
                    '/jointset/back/lumbar_rotation/value',
                    '/jointset/subtalar_l/subtalar_angle_l/value',
                    '/jointset/subtalar_r/subtalar_angle_r/value',
                    '/jointset/mtp_l/mtp_angle_l/value',
                    '/jointset/mtp_r/mtp_angle_r/value']
    
    #Loop through and remove values
    for rr in range(0,len(removeLabels)):
        origStorage.removeColumn(removeLabels[rr])
        
    #Get current column labels 
    angleNames = list(origStorage.getColumnLabels())
    
    # #Rename walker knee to standard knee joint
    # angleNames = [sub.replace('walker_knee_', 'knee_') for sub in angleNames] 
    
    #Output as temp .sto to get back as storage
    stoAdapter = osim.STOFileAdapter()
    stoAdapter.write(origStorage,'temp.sto')
    newStorage = osim.Storage('temp.sto')
    
    #Rename column labels in storage
    #Create array string
    colLabels = osim.ArrayStr()
    #Append time label
    colLabels.append('time')
    #Append the state labels
    for aa in range(0,len(angleNames)):
        colLabels.append(angleNames[aa])
    #Rename labels in storage
    newStorage.setColumnLabels(colLabels)
    
    # #Figure out the state index of the knee angles
    # rKneeInd = angleNames.index('/jointset/knee_r/knee_angle_r/value')
    # lKneeInd = angleNames.index('/jointset/knee_l/knee_angle_l/value')

    # #Invert knee angles to match different model
    # #Loop through the length of the data
    # for cc in range(0,newStorage.getSize()):
    #     #Get the current knee angle values
    #     rKneeVal = newStorage.getStateVector(cc).getData().getitem(rKneeInd)
    #     lKneeVal = newStorage.getStateVector(cc).getData().getitem(lKneeInd)
    #     #Replace with inverted value
    #     newStorage.getStateVector(cc).getData().setitem(rKneeInd,rKneeVal*-1)
    #     newStorage.getStateVector(cc).getData().setitem(lKneeInd,lKneeVal*-1)
    
    #Cleanup temp.sto
    os.remove('temp.sto')

    #Write the states storage object to file
    newStorage.printToXML(outputFileName)

# %% addCoordinateActuator

def addCoordinateActuator(osimModel = None,
                          coordName = None,
                          optForce = 1000,
                          controlVals = None):
    
    # Convenience function for adding torque coordinate actuators to model
    #
    # Input:    osimModel - model object to add actuator to
    #           coordName - name of coordinate to actuate
    #           optForce - optimal force value for coordinate actuator
    #           controlVals - min and maximum control values for actuator (if None Inf is used)
    
    #Input checks
    if osimModel is None:
        raise ValueError('A model object is required.')
        
    if coordName is None:
        raise ValueError('A coordinate name is required.')
    
    #Get coordinate set from model
    coordSet = osimModel.updCoordinateSet()
    
    #Create actuator
    actu = osim.CoordinateActuator()
    
    #Set actuator name
    actu.setName('tau_'+coordName)
    
    #Set coordinate for actuator
    actu.setCoordinate(coordSet.get(coordName))
    
    #Set optimal force for actuator
    actu.setOptimalForce(optForce)
    
    #Set max and min controls
    if controlVals is None:
        actu.setMinControl(np.inf*-1)
        actu.setMaxControl(np.inf)
    else:
        actu.setMinControl(controlVals[0])
        actu.setMaxControl(controlVals[1])
        
    #Add actuator to model
    osimModel.updForceSet().cloneAndAppend(actu)
    
# %% mm_to_m

def mm_to_m(table, label):
    
    #Function for converting mm to m units
    
    c = table.updDependentColumn(label)
    for i in range(c.size()):
        c[i] = osim.Vec3(c[i][0] * 0.001, c[i][1] * 0.001, c[i][2] * 0.001)

# %% readSTO

def readSTO(fileName):
    
    #Read storage
    sto = osim.Storage(fileName)
        
    #Get column list
    labels = osimArrayToList(sto.getColumnLabels())

    #Get time
    time = osim.ArrayDouble()
    sto.getTimeColumn(time)  
    time = osimArrayToList(time)
    
    #Get data
    data = []
    for i in range(sto.getSize()):
        temp = osimArrayToList(sto.getStateVector(i).getData())
        temp.insert(0, time[i])
        data.append(temp)

    df = pd.DataFrame(data, columns=labels)
    df.index = df.time
    
    return df

# %% osimArrayToList

def osimArrayToList(array):
    temp = []
    for i in range(array.getSize()):
        temp.append(array.get(i))

    return temp

# %% inverseSolutionToTrackGuess
    
def inverseSolutionToTrackGuess(guessFile = None, mocoSolver = None):
        
    # Convenience function for fixing a guess file that is generated from an inverse
    # solution to be placed in a track solver. There appears to be a slight mismatch
    # between the guess from a MocoTrack to the solution from a MocoInverse. This
    # can be fixed by manually placing the parameters of the solution into a randomly
    # created guess.
    #
    # Input:    guessFile - string to path of guess file
    #           mocoSolver - solver object to place guess in
    # 
    # Note: this function will only work if the guess file matches the inputs
    # to the solver exactly. If they don't match things won't work well...
    
    ##### TODO: this currently works for baseline sims, but may need some options if re-using in other processes...
    
    #Check for appropriate inputs
    if guessFile is None or mocoSolver is None:
        raise ValueError('A guess file and linked Moco Solver are needed in fixGuessFile')
        
    #Grab the file as a Moco tracjectory
    mocoTraj = osim.MocoTrajectory(guessFile)

    #Create random guess using the current solver
    randTraj = mocoSolver.createGuess()
    
    #Resample the randomly created guess if it doesn't match the guess file
    if randTraj.getNumTimes() != mocoTraj.getNumTimes():
        randTraj.resampleWithNumTimes(mocoTraj.getNumTimes())
    
    #Set time for random guess to original trajectory
    randTraj.setTime(mocoTraj.getTime())

    #Refresh states in random guess from the moco trajectory
    randTraj.setStatesTrajectory(mocoTraj.exportToStatesTable())

    #Refresh controls from moco trajectory
    controlNames = list(mocoTraj.getControlNames())
    #Loop through and set the controls in the guess
    for cc in range(0,len(controlNames)):
        randTraj.setControl(controlNames[cc],mocoTraj.getControlMat(controlNames[cc]))

    #Refresh multipliers from moco trajectory
    multiplierNames = list(mocoTraj.getMultiplierNames())
    for cc in range(0,len(multiplierNames)):
        randTraj.setMultiplier(multiplierNames[cc],mocoTraj.getMultiplierMat(multiplierNames[cc]))
        
    #Set the guess in the input solver
    mocoSolver.setGuess(randTraj)   
    
    
