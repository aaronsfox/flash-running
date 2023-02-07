# -*- coding: utf-8 -*-
"""

@author: 
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    Build of functions to assist with running simulations
    
"""


import opensim as osim
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.signal import butter, filtfilt
import os
import shutil

# %% Function to add set of torque actuators to model

def addTorqueActuators(osimModel = None,
                       optForces = None):
    
    """
    
    Convenience function for adding series of torque actuators to model
    
    Input:    osimModel - OpenSim model object for use
              optForces - dict of coordinates and their associated optimal forces to add
              
    Output:   osimModel - updated torque driven model
                  
    """
    
    #Check inputs
    if osimModel is None or optForces is None:
            raise ValueError('All inputs for this function are required!')
    
    #Remove the original lumbar actuators to not apply force
    #Get the force indices to remove from the model
    forceRemove = []
    for forceInd in range(osimModel.updForceSet().getSize()):
        if 'lumbar' in osimModel.updForceSet().get(forceInd).getName():
            forceRemove.append(forceInd)
    #Each time the force is removed the indices reduce, so need to account for this
    for ind in range(len(forceRemove)):
        osimModel.updForceSet().remove(forceRemove[ind]-ind)
    
    #Intialise model system
    osimModel.initSystem()
    
    #Get coordinate list
    coordinatesList = list(optForces.keys())
    
    #Get coordinate set
    coordSet = osimModel.getCoordinateSet()
    
    #Loop through coordinates and add actuators
    for coordinate in coordinatesList:
        #Create actuator
        actu = osim.CoordinateActuator()
        #Set name
        actu.setName(f'{coordinate}_actuator')
        #Set coordinate
        actu.setCoordinate(coordSet.get(coordinate))
        #Set optimal force
        actu.setOptimalForce(optForces[coordinate])
        #Set min and max control
        actu.setMinControl(-1)
        actu.setMaxControl(1)
        #Append to model force set
        osimModel.updForceSet().cloneAndAppend(actu)
        # #Append to model components
        # osimModel.addComponent(actu)
    
    #Finalise model connections
    osimModel.finalizeConnections()
    
    #Return model
    return osimModel

# %% ----- End of osimFunctions.py -----