# -*- coding: utf-8 -*-
"""

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This script tests running a simulation of sprint running using OpenSim AD as
    outlined in the TGCS 2023 OpenSim workshop. See comments throughout code for
    further info.


"""

# %% Import packages

import os
import sys
import logging
import opensim as osim
import numpy as np
import platform
import shutil
import importlib
import casadi as ca
import pandas as pd

# %% Define some basic functions

#Storage file to dataframe.
def storage_to_dataframe(storage_file, headers):
    # Extract data
    data = storage_to_numpy(storage_file)
    out = pd.DataFrame(data=data['time'], columns=['time'])    
    for count, header in enumerate(headers):
        out.insert(count + 1, header, data[header])    
    
    return out

#Storage file to numpy array.
def storage_to_numpy(storage_file, excess_header_entries=0):
    """Returns the data from a storage file in a numpy format. Skips all lines
    up to and including the line that says 'endheader'.
    Parameters
    ----------
    storage_file : str
        Path to an OpenSim Storage (.sto) file.
    Returns
    -------
    data : np.ndarray (or numpy structure array or something?)
        Contains all columns from the storage file, indexable by column name.
    excess_header_entries : int, optional
        If the header row has more names in it than there are data columns.
        We'll ignore this many header row entries from the end of the header
        row. This argument allows for a hacky fix to an issue that arises from
        Static Optimization '.sto' outputs.
    Examples
    --------
    Columns from the storage file can be obtained as follows:
        >>> data = storage2numpy('<filename>')
        >>> data['ground_force_vy']
    """
    # What's the line number of the line containing 'endheader'?
    f = open(storage_file, 'r')

    header_line = False
    for i, line in enumerate(f):
        if header_line:
            column_names = line.split()
            break
        if line.count('endheader') != 0:
            line_number_of_line_containing_endheader = i + 1
            header_line = True
    f.close()

    # With this information, go get the data.
    if excess_header_entries == 0:
        names = True
        skip_header = line_number_of_line_containing_endheader
    else:
        names = column_names[:-excess_header_entries]
        skip_header = line_number_of_line_containing_endheader + 1
    data = np.genfromtxt(storage_file, names=names,
            skip_header=skip_header)

    return data

#Settings
def get_setup(motion_type):

    setups = {}   
    setups['other'] = {
        'ipopt_tolerance': 3,
        'weights': {
            'positionTrackingTerm': 100,
            'velocityTrackingTerm': 10,
            'accelerationTrackingTerm': 50,
            'activationTerm': 10,
            'armExcitationTerm': 0.001,
            'lumbarExcitationTerm': 0.001,
            'jointAccelerationTerm': 0.001,
            'activationDtTerm': 0.001,
            'forceDtTerm': 0.001},            
        'coordinates_toTrack': {
            'pelvis_tilt': {"weight": 10},
            'pelvis_list': {"weight": 10},
            'pelvis_rotation': {"weight": 10},
            'pelvis_tx': {"weight": 10},
            'pelvis_ty': {"weight": 10},
            'pelvis_tz': {"weight": 10}, 
            'hip_flexion_l': {"weight": 20},
            'hip_adduction_l': {"weight": 10},
            'hip_rotation_l': {"weight": 1},
            'hip_flexion_r': {"weight": 20},
            'hip_adduction_r': {"weight": 10},
            'hip_rotation_r': {"weight": 1},
            'knee_angle_l': {"weight": 10},
            'knee_angle_r': {"weight": 10},
            'ankle_angle_l': {"weight": 10},
            'ankle_angle_r': {"weight": 10},
            'subtalar_angle_l': {"weight": 10},
            'subtalar_angle_r': {"weight": 10},
            'lumbar_extension': {"weight": 10},
            'lumbar_bending': {"weight": 10},
            'lumbar_rotation': {"weight": 10},
            'arm_flex_l': {"weight": 10},
            'arm_add_l': {"weight": 10},
            'arm_rot_l': {"weight": 10},
            'arm_flex_r': {"weight": 10},
            'arm_add_r': {"weight": 10},
            'arm_rot_r': {"weight": 10},
            'elbow_flex_l': {"weight": 10},
            'elbow_flex_r': {"weight": 10},
            'pro_sup_l': {"weight": 10},
            'pro_sup_r': {"weight": 10}},
        'coordinate_constraints': {
            'pelvis_tx': {"env_bound": 0.1}},
        'ignorePassiveFiberForce': True,
        'filter_Qs_toTrack': True,
        'cutoff_freq_Qs': 30,
        'filter_Qds_toTrack': True,
        'cutoff_freq_Qds': 30,
        'filter_Qdds_toTrack': True,
        'cutoff_freq_Qdds': 30,
        'splineQds': True,
        'meshDensity': 100,
        'yCalcnToes': True}
    
    setups['running'] = {
        'ipopt_tolerance': 3,
        'weights': {
            'positionTrackingTerm': 100,
            'velocityTrackingTerm': 10,
            'accelerationTrackingTerm': 50,
            'activationTerm': 10,
            'armExcitationTerm': 0.001,
            'lumbarExcitationTerm': 0.001,
            'jointAccelerationTerm': 0.001,
            'activationDtTerm': 0.001,
            'forceDtTerm': 0.001},            
        'coordinates_toTrack': {
            'pelvis_tilt': {"weight": 10},
            'pelvis_list': {"weight": 10},
            'pelvis_rotation': {"weight": 10},
            'pelvis_tx': {"weight": 10},
            'pelvis_ty': {"weight": 10},
            'pelvis_tz': {"weight": 10}, 
            'hip_flexion_l': {"weight": 20},
            'hip_adduction_l': {"weight": 10},
            'hip_rotation_l': {"weight": 1},
            'hip_flexion_r': {"weight": 20},
            'hip_adduction_r': {"weight": 10},
            'hip_rotation_r': {"weight": 1},
            'knee_angle_l': {"weight": 10},
            'knee_angle_r': {"weight": 10},
            'ankle_angle_l': {"weight": 10},
            'ankle_angle_r': {"weight": 10},
            'subtalar_angle_l': {"weight": 10},
            'subtalar_angle_r': {"weight": 10},
            'lumbar_extension': {"weight": 10},
            'lumbar_bending': {"weight": 10},
            'lumbar_rotation': {"weight": 10},
            'arm_flex_l': {"weight": 10},
            'arm_add_l': {"weight": 10},
            'arm_rot_l': {"weight": 10},
            'arm_flex_r': {"weight": 10},
            'arm_add_r': {"weight": 10},
            'arm_rot_r': {"weight": 10},
            'elbow_flex_l': {"weight": 10},
            'elbow_flex_r': {"weight": 10},
            'pro_sup_l': {"weight": 10},
            'pro_sup_r': {"weight": 10}},
        'coordinate_constraints': {
            'pelvis_tx': {"env_bound": 0.1}},
        'ignorePassiveFiberForce': True,
        'filter_Qs_toTrack': True,
        'cutoff_freq_Qs': 12,
        'filter_Qds_toTrack': True,
        'cutoff_freq_Qds': 12,
        'filter_Qdds_toTrack': True,
        'cutoff_freq_Qdds': 12,
        'splineQds': True,
        'meshDensity': 100,
        'yCalcnToes': True}
    
    setups['running_torque_driven'] = {
        'ipopt_tolerance': 3,
        'weights': {
            'positionTrackingTerm': 100,
            'velocityTrackingTerm': 10,
            'accelerationTrackingTerm': 50,
            'armExcitationTerm': 0.001,
            'lumbarExcitationTerm': 0.001,
            'jointAccelerationTerm': 0.001,
            'coordinateExcitationTerm': 10},            
        'coordinates_toTrack': {
            'pelvis_tilt': {"weight": 10},
            'pelvis_list': {"weight": 10},
            'pelvis_rotation': {"weight": 10},
            'pelvis_tx': {"weight": 10},
            'pelvis_ty': {"weight": 10},
            'pelvis_tz': {"weight": 10}, 
            'hip_flexion_l': {"weight": 20},
            'hip_adduction_l': {"weight": 10},
            'hip_rotation_l': {"weight": 1},
            'hip_flexion_r': {"weight": 20},
            'hip_adduction_r': {"weight": 10},
            'hip_rotation_r': {"weight": 1},
            'knee_angle_l': {"weight": 10},
            'knee_angle_r': {"weight": 10},
            'ankle_angle_l': {"weight": 10},
            'ankle_angle_r': {"weight": 10},
            'subtalar_angle_l': {"weight": 10},
            'subtalar_angle_r': {"weight": 10},
            'lumbar_extension': {"weight": 10},
            'lumbar_bending': {"weight": 10},
            'lumbar_rotation': {"weight": 10},
            'arm_flex_l': {"weight": 10},
            'arm_add_l': {"weight": 10},
            'arm_rot_l': {"weight": 10},
            'arm_flex_r': {"weight": 10},
            'arm_add_r': {"weight": 10},
            'arm_rot_r': {"weight": 10},
            'elbow_flex_l': {"weight": 10},
            'elbow_flex_r': {"weight": 10},
            'pro_sup_l': {"weight": 10},
            'pro_sup_r': {"weight": 10}},
        'coordinate_constraints': {
            'pelvis_tx': {"env_bound": 0.1}},
        'ignorePassiveFiberForce': True,
        'filter_Qs_toTrack': True,
        'cutoff_freq_Qs': 12,
        'filter_Qds_toTrack': True,
        'cutoff_freq_Qds': 12,
        'filter_Qdds_toTrack': True,
        'cutoff_freq_Qdds': 12,
        'splineQds': True,
        'meshDensity': 100,
        'yCalcnToes': True,
        'torque_driven_model': True,
        'coordinate_optimal_forces': {
            'hip_flexion_r': 400,
            'hip_flexion_l': 400}}
    
    setups['walking'] = {
        'ipopt_tolerance': 3,
        'weights': {
            'positionTrackingTerm': 10,
            'velocityTrackingTerm': 10,
            'accelerationTrackingTerm': 50,
            'activationTerm': 1,
            'armExcitationTerm': 0.001,
            'lumbarExcitationTerm': 0.001,
            'jointAccelerationTerm': 0.001,
            'activationDtTerm': 0.001,
            'forceDtTerm': 0.001},            
        'coordinates_toTrack': {
            'pelvis_tilt': {"weight": 10},
            'pelvis_list': {"weight": 1},
            'pelvis_rotation': {"weight": 1},
            'pelvis_tx': {"weight": 1},
            'pelvis_ty': {"weight": 1},
            'pelvis_tz': {"weight": 1}, 
            'hip_flexion_l': {"weight": 10},
            'hip_adduction_l': {"weight": 1},
            'hip_rotation_l': {"weight": 1},
            'hip_flexion_r': {"weight": 10},
            'hip_adduction_r': {"weight": 1},
            'hip_rotation_r': {"weight": 1},
            'knee_angle_l': {"weight": 10},
            'knee_angle_r': {"weight": 10},
            'ankle_angle_l': {"weight": 10},
            'ankle_angle_r': {"weight": 10},
            'subtalar_angle_l': {"weight": 1},
            'subtalar_angle_r': {"weight": 1},
            'lumbar_extension': {"weight": 10},
            'lumbar_bending': {"weight": 1},
            'lumbar_rotation': {"weight": 1},
            'arm_flex_l': {"weight": 1},
            'arm_add_l': {"weight": 1},
            'arm_rot_l': {"weight": 1},
            'arm_flex_r': {"weight": 1},
            'arm_add_r': {"weight": 1},
            'arm_rot_r': {"weight": 1},
            'elbow_flex_l': {"weight": 1},
            'elbow_flex_r': {"weight": 1},
            'pro_sup_l': {"weight": 1},
            'pro_sup_r': {"weight": 1}},            
        'coordinate_constraints': {
            'pelvis_ty': {"env_bound": 0.1},
            'pelvis_tx': {"env_bound": 0.1}},
        'enableLimitTorques': True,
        'filter_Qs_toTrack': True,
        'cutoff_freq_Qs': 6,
        'meshDensity': 100}
    
    setups['drop_jump'] = {
        'weights': {
            'positionTrackingTerm': 50,
            'velocityTrackingTerm': 10,
            'accelerationTrackingTerm': 50,
            'activationTerm': 1,
            'armExcitationTerm': 0.001,
            'lumbarExcitationTerm': 0.001,
            'jointAccelerationTerm': 0.001,
            'activationDtTerm': 0.001,
            'forceDtTerm': 0.001},            
        'coordinates_toTrack': {
            'pelvis_tilt': {"weight": 10},
            'pelvis_list': {"weight": 1},
            'pelvis_rotation': {"weight": 1},
            'pelvis_tx': {"weight": 1},
            'pelvis_ty': {"weight": 10},
            'pelvis_tz': {"weight": 1}, 
            'hip_flexion_l': {"weight": 10},
            'hip_adduction_l': {"weight": 1},
            'hip_rotation_l': {"weight": 1},
            'hip_flexion_r': {"weight": 10},
            'hip_adduction_r': {"weight": 1},
            'hip_rotation_r': {"weight": 1},
            'knee_angle_l': {"weight": 10},
            'knee_angle_r': {"weight": 10},
            'ankle_angle_l': {"weight": 10},
            'ankle_angle_r': {"weight": 10},
            'subtalar_angle_l': {"weight": 1},
            'subtalar_angle_r': {"weight": 1},
            'lumbar_extension': {"weight": 10},
            'lumbar_bending': {"weight": 1},
            'lumbar_rotation': {"weight": 1},
            'arm_flex_l': {"weight": 50},
            'arm_add_l': {"weight": 50},
            'arm_rot_l': {"weight": 50},
            'arm_flex_r': {"weight": 50},
            'arm_add_r': {"weight": 50},
            'arm_rot_r': {"weight": 50},
            'elbow_flex_l': {"weight": 50},
            'elbow_flex_r': {"weight": 50},
            'pro_sup_l': {"weight": 50},
            'pro_sup_r': {"weight": 50}},            
        'coordinate_constraints': {
            'pelvis_ty': {"env_bound": 0.02},
            'pelvis_tx': {"env_bound": 0.02},
            'pelvis_tz': {"env_bound": 0.02}},
        'ignorePassiveFiberForce': True,
        'filter_Qs_toTrack': True,
        'cutoff_freq_Qs': 30,
        'filter_Qds_toTrack': True,
        'cutoff_freq_Qds': 30,
        'filter_Qdds_toTrack': True,
        'cutoff_freq_Qdds': 30,
        'splineQds': True,
        'meshDensity': 100}
    
    setups['sit_to_stand'] = {
        'ipopt_tolerance': 3,
        'weights': {
            'positionTrackingTerm': 50,
            'velocityTrackingTerm': 10,
            'accelerationTrackingTerm': 50,
            'activationTerm': 100,
            'armExcitationTerm': 0.001,
            'lumbarExcitationTerm': 0.001,
            'jointAccelerationTerm': 0.001,
            'activationDtTerm': 0.001,
            'forceDtTerm': 0.001,
            'reserveActuatorTerm': 0.001,
            },            
        'coordinates_toTrack': {
            'pelvis_tilt': {"weight": 100},
            'pelvis_list': {"weight": 10},
            'pelvis_rotation': {"weight": 1},
            'pelvis_tx': {"weight": 100},
            'pelvis_ty': {"weight": 10},
            'pelvis_tz': {"weight": 100}, 
            'hip_flexion_l': {"weight": 100},
            'hip_adduction_l': {"weight": 20},
            'hip_rotation_l': {"weight": 1},
            'hip_flexion_r': {"weight": 100},
            'hip_adduction_r': {"weight": 20},
            'hip_rotation_r': {"weight": 1},
            'knee_angle_l': {"weight": 100},
            'knee_angle_r': {"weight": 100},
            'ankle_angle_l': {"weight": 100},
            'ankle_angle_r': {"weight": 100},
            'subtalar_angle_l': {"weight": 20},
            'subtalar_angle_r': {"weight": 20},
            'lumbar_extension': {"weight": 100},
            'lumbar_bending': {"weight": 20},
            'lumbar_rotation': {"weight": 20},
            'arm_flex_l': {"weight": 50},
            'arm_add_l': {"weight": 10},
            'arm_rot_l': {"weight": 10},
            'arm_flex_r': {"weight": 50},
            'arm_add_r': {"weight": 10},
            'arm_rot_r': {"weight": 10},
            'elbow_flex_l': {"weight": 10},
            'elbow_flex_r': {"weight": 10},
            'pro_sup_l': {"weight": 10},
            'pro_sup_r': {"weight": 10}},            
        'coordinate_constraints': {
            'pelvis_ty': {"env_bound": 0.1},
            'pelvis_tx': {"env_bound": 0.1}},       
        'withReserveActuators': True,
        'reserveActuatorCoordinates': {
            'hip_rotation_l': 30, 'hip_rotation_r': 30},
        'periodicConstraints': {'coordinateValues': ['lowerLimbJoints']},
        'ignorePassiveFiberForce': True,
        'filter_Qs_toTrack': True,
        'cutoff_freq_Qs': 4,
        'filter_Qds_toTrack': True,
        'cutoff_freq_Qds': 4,
        'filter_Qdds_toTrack': True,
        'cutoff_freq_Qdds': 4,
        'splineQds': True,
        'meshDensity': 50}
    
    setups['squats'] = {
        'ipopt_tolerance': 3,
        'weights': {
            'positionTrackingTerm': 50,
            'velocityTrackingTerm': 10,
            'accelerationTrackingTerm': 50,
            'activationTerm': 100,
            'armExcitationTerm': 0.001,
            'lumbarExcitationTerm': 0.001,
            'jointAccelerationTerm': 0.001,
            'activationDtTerm': 0.001,
            'forceDtTerm': 0.001,
            'reserveActuatorTerm': 0.001},            
        'coordinates_toTrack': {
            'pelvis_tilt': {"weight": 100},
            'pelvis_list': {"weight": 10},
            'pelvis_rotation': {"weight": 1},
            'pelvis_tx': {"weight": 100},
            'pelvis_ty': {"weight": 10},
            'pelvis_tz': {"weight": 100}, 
            'hip_flexion_l': {"weight": 100},
            'hip_adduction_l': {"weight": 20},
            'hip_rotation_l': {"weight": 1},
            'hip_flexion_r': {"weight": 100},
            'hip_adduction_r': {"weight": 20},
            'hip_rotation_r': {"weight": 1},
            'knee_angle_l': {"weight": 100},
            'knee_angle_r': {"weight": 100},
            'ankle_angle_l': {"weight": 100},
            'ankle_angle_r': {"weight": 100},
            'subtalar_angle_l': {"weight": 20},
            'subtalar_angle_r': {"weight": 20},
            'lumbar_extension': {"weight": 100},
            'lumbar_bending': {"weight": 20},
            'lumbar_rotation': {"weight": 20},
            'arm_flex_l': {"weight": 50},
            'arm_add_l': {"weight": 10},
            'arm_rot_l': {"weight": 10},
            'arm_flex_r': {"weight": 50},
            'arm_add_r': {"weight": 10},
            'arm_rot_r': {"weight": 10},
            'elbow_flex_l': {"weight": 10},
            'elbow_flex_r': {"weight": 10},
            'pro_sup_l': {"weight": 10},
            'pro_sup_r': {"weight": 10}},            
        'coordinate_constraints': {
            'pelvis_ty': {"env_bound": 0.1},
            'pelvis_tx': {"env_bound": 0.1}},
        'withReserveActuators': True,
        'reserveActuatorCoordinates': {
            'hip_rotation_l': 30, 'hip_rotation_r': 30},
        'periodicConstraints': {'coordinateValues': ['lowerLimbJoints'],
                                'coordinateSpeeds': ['lowerLimbJoints'],
                                'muscleActivationsForces': ['all'],
                                'lumbarJointActivations': ['all']},
        'ignorePassiveFiberForce': True,
        'filter_Qs_toTrack': True,
        'cutoff_freq_Qs': 4,
        'filter_Qds_toTrack': True,
        'cutoff_freq_Qds': 4,
        'filter_Qdds_toTrack': True,
        'cutoff_freq_Qdds': 4,
        'splineQds': True,
        'heel_vGRF_threshold': 5,
        'meshDensity': 50}
        
    setups['jumping'] = {
        'weights': {
            'positionTrackingTerm': 100,
            'velocityTrackingTerm': 10,
            'accelerationTrackingTerm': 50,
            'activationTerm': 1,
            'armExcitationTerm': 0.001,
            'lumbarExcitationTerm': 0.001,
            'jointAccelerationTerm': 0.001,
            'activationDtTerm': 0.001,
            'forceDtTerm': 0.001},            
        'coordinates_toTrack': {
            'pelvis_tilt': {"weight": 10},
            'pelvis_list': {"weight": 10},
            'pelvis_rotation': {"weight": 10},
            'pelvis_tx': {"weight": 10},
            'pelvis_ty': {"weight": 100},
            'pelvis_tz': {"weight": 10}, 
            'hip_flexion_l': {"weight": 20},
            'hip_adduction_l': {"weight": 10},
            'hip_rotation_l': {"weight": 10},
            'hip_flexion_r': {"weight": 20},
            'hip_adduction_r': {"weight": 10},
            'hip_rotation_r': {"weight": 10},
            'knee_angle_l': {"weight": 10},
            'knee_angle_r': {"weight": 10},
            'ankle_angle_l': {"weight": 10},
            'ankle_angle_r': {"weight": 10},
            'subtalar_angle_l': {"weight": 10},
            'subtalar_angle_r': {"weight": 10},
            'lumbar_extension': {"weight": 10},
            'lumbar_bending': {"weight": 10},
            'lumbar_rotation': {"weight": 10},
            'arm_flex_l': {"weight": 100},
            'arm_add_l': {"weight": 100},
            'arm_rot_l': {"weight": 100},
            'arm_flex_r': {"weight": 100},
            'arm_add_r': {"weight": 100},
            'arm_rot_r': {"weight": 100},
            'elbow_flex_l': {"weight": 100},
            'elbow_flex_r': {"weight": 100},
            'pro_sup_l': {"weight": 100},
            'pro_sup_r': {"weight": 100}},
        'coordinate_constraints': {
            'pelvis_tx': {"env_bound": 0.1},
            'pelvis_ty': {"env_bound": 0.1}},
        'ignorePassiveFiberForce': True,
        'filter_Qs_toTrack': True,
        'cutoff_freq_Qs': 10,
        'filter_Qds_toTrack': True,
        'cutoff_freq_Qds': 10,
        'filter_Qdds_toTrack': True,
        'cutoff_freq_Qdds': 10,
        'splineQds': True,
        'meshDensity': 50,
        'yCalcnToes': True,
        }
    
    setups['my_periodic_running'] = {
        'ipopt_tolerance': 3,
        'weights': {
            'positionTrackingTerm': 100,
            'velocityTrackingTerm': 10,
            'accelerationTrackingTerm': 50,
            'activationTerm': 10,
            'armExcitationTerm': 0.001,
            'lumbarExcitationTerm': 0.001,
            'jointAccelerationTerm': 0.001,
            'activationDtTerm': 0.001,
            'forceDtTerm': 0.001},            
        'coordinates_toTrack': {
            'pelvis_tilt': {"weight": 10},
            'pelvis_list': {"weight": 10},
            'pelvis_rotation': {"weight": 10},
            'pelvis_tx': {"weight": 10},
            'pelvis_ty': {"weight": 10},
            'pelvis_tz': {"weight": 10}, 
            'hip_flexion_l': {"weight": 20},
            'hip_adduction_l': {"weight": 10},
            'hip_rotation_l': {"weight": 1},
            'hip_flexion_r': {"weight": 20},
            'hip_adduction_r': {"weight": 10},
            'hip_rotation_r': {"weight": 1},
            'knee_angle_l': {"weight": 10},
            'knee_angle_r': {"weight": 10},
            'ankle_angle_l': {"weight": 10},
            'ankle_angle_r': {"weight": 10},
            'subtalar_angle_l': {"weight": 10},
            'subtalar_angle_r': {"weight": 10},
            'lumbar_extension': {"weight": 10},
            'lumbar_bending': {"weight": 10},
            'lumbar_rotation': {"weight": 10},
            'arm_flex_l': {"weight": 10},
            'arm_add_l': {"weight": 10},
            'arm_rot_l': {"weight": 10},
            'arm_flex_r': {"weight": 10},
            'arm_add_r': {"weight": 10},
            'arm_rot_r': {"weight": 10},
            'elbow_flex_l': {"weight": 10},
            'elbow_flex_r': {"weight": 10},
            'pro_sup_l': {"weight": 10},
            'pro_sup_r': {"weight": 10}},
        'coordinate_constraints': {
            'pelvis_tx': {"env_bound": 0.1}},
        'periodicConstraints': {
            # All lower limb coordinates but pelvis_tx.
            'coordinateValues': ['pelvis_tilt', 'pelvis_list', 'pelvis_rotation', 
                   'pelvis_ty', 'pelvis_tz', 'hip_flexion_l', 
                   'hip_adduction_l', 'hip_rotation_l', 'hip_flexion_r',
                   'hip_adduction_r', 'hip_rotation_r', 'knee_angle_l',
                   'knee_angle_r', 'ankle_angle_l', 'ankle_angle_r', 
                   'subtalar_angle_l', 'subtalar_angle_r', 'mtp_angle_l',
                   'mtp_angle_r', 'lumbar_extension', 'lumbar_bending',
                   'lumbar_rotation'],
            'coordinateSpeeds': ['lowerLimbJoints'],
            'muscleActivationsForces': ['all'],
            'lumbarJointActivations': ['all']},
        'ignorePassiveFiberForce': True,
        'filter_Qs_toTrack': True,
        'cutoff_freq_Qs': 12,
        'filter_Qds_toTrack': True,
        'cutoff_freq_Qds': 12,
        'filter_Qdds_toTrack': True,
        'cutoff_freq_Qdds': 12,
        'splineQds': True,
        'meshDensity': 100,
        'yCalcnToes': True}

    return setups[motion_type]

# %% Set-up

# #Set matplotlib parameters
# from matplotlib import rcParams
# # rcParams['font.family'] = 'sans-serif'
# rcParams['font.sans-serif'] = 'Arial'
# rcParams['font.weight'] = 'bold'
# rcParams['axes.labelsize'] = 12
# rcParams['axes.titlesize'] = 16
# rcParams['axes.linewidth'] = 1.5
# rcParams['axes.labelweight'] = 'bold'
# rcParams['legend.fontsize'] = 10
# rcParams['xtick.major.width'] = 1.5
# rcParams['ytick.major.width'] = 1.5
# rcParams['legend.framealpha'] = 0.0
# rcParams['savefig.dpi'] = 300
# rcParams['savefig.format'] = 'pdf'

# Insert the name of the trial you want to simulate.
trial_name = 'coordinates'

# Insert the type of activity you want to simulate. We have pre-defined settings
# for different activities (more details above). Visit 
# ./UtilsDynamicSimulations/OpenSimAD/settingsOpenSimAD.py to see all available
# activities and their settings. If your activity is not in the list, select
# 'other' to use generic settings or set your own settings.
motion_type = 'running_torque_driven' ###note this was changed

# Insert the time interval you want to simulate. It is recommended to simulate
# trials shorter than 2s (more details above). Set to [] to simulate full trial.
# Note the time window is extracted from some earlier simulation tests
time_window = [0.352, 0.582]

# Insert the speed of the treadmill in m/s. A positive value indicates that the
# subject is moving forward. You should ignore this parameter or set it to 0 if
# the trial was not measured on a treadmill.
treadmill_speed = 0

# %% Process inputs for OpenSim AD

"""

The following sections of code manually replicates the processInputsOpenSimAD
function included with the example script in the opencap-processing repository.
The processInputsOpenSimAD functions is in the utilsOpenSimAD script also included.
This therefore could once again be wrapped up in processing at a later time.
Included in these sections:
    > Scale model
    > Adjusts wrapping surfaces
    > Adds foot ground contacts
    > Generates external functions for OpenSim AD

"""

# %% Scale the JA1 model to the LaiUhlrich model

#Set the XYZ marker pairs for each body in the model
scaleMarkerPairs = {
    'pelvis': [
        [['RASI', 'SACR'], ['LASI', 'SACR']],
        [['RASI', 'LASI']],
        [['RASI', 'LASI']]
        ],
    'femur_r': [
        [['RLEPI', 'RMEPI']],
        [['RASI', 'RLEPI']],
        [['RLEPI', 'RMEPI']]
        ],
    'tibia_r': [
        [['RLMAL', 'RMMAL']],
        [['RLEPI', 'RLMAL'], ['RMEPI', 'RMMAL']],
        [['RLEPI', 'RMEPI']]
        ],
    'patella_r': [
        [['RLEPI', 'RMEPI']],
        [['RLEPI', 'RMEPI']],
        [['RLEPI', 'RMEPI']]
        ],
    'talus_r': [
        [['RMMAL', 'RLMAL']],
        [['RMMAL', 'RLMAL']],
        [['RMMAL', 'RLMAL']],
        ],
    'calcn_r': [
        [['RHEEL', 'RLMAL']],
        [['RHEEL', 'RLMAL']],
        [['RHEEL', 'RLMAL']],
        ],
    'toes_r': [
        [['RHEEL', 'RTOE']],
        [['RHEEL', 'RTOE']],
        [['RP1MT', 'RP5MT']],
        ],
    'femur_l': [
        [['LLEPI', 'LMEPI']],
        [['LASI', 'LLEPI']],
        [['LLEPI', 'LMEPI']]
        ],
    'tibia_l': [
        [['LLMAL', 'LMMAL']],
        [['LLEPI', 'LLMAL'], ['LMEPI', 'LMMAL']],
        [['LLEPI', 'LMEPI']]
        ],
    'patella_l': [
        [['LLEPI', 'LMEPI']],
        [['LLEPI', 'LMEPI']],
        [['LLEPI', 'LMEPI']]
        ],
    'talus_l': [
        [['LMMAL', 'LLMAL']],
        [['LMMAL', 'LLMAL']],
        [['LMMAL', 'LLMAL']],
        ],
    'calcn_l': [
        [['LHEEL', 'LLMAL']],
        [['LHEEL', 'LLMAL']],
        [['LHEEL', 'LLMAL']],
        ],
    'toes_l': [
        [['LHEEL', 'LTOE']],
        [['LHEEL', 'LTOE']],
        [['LP1MT', 'LP5MT']],
        ],
    'torso': [
        [['MAN', 'C7']],
        [['RSH', 'RASI'],['LSH', 'LASI']],
        [['RSH', 'LSH']],
        ],
    'humerus_r': [
        [['RSH', 'RELB']],
        [['RSH', 'RELB']],
        [['RSH', 'RELB']],
        ],
    'ulna_r': [
        [['RELB', 'RWR']],
        [['RELB', 'RWR']],
        [['RELB', 'RWR']],
        ],
    'radius_r': [
        [['RELB', 'RWR']],
        [['RELB', 'RWR']],
        [['RELB', 'RWR']],
        ],
    'humerus_l': [
        [['LSH', 'LELB']],
        [['LSH', 'LELB']],
        [['LSH', 'LELB']],
        ],
    'ulna_l': [
        [['LELB', 'LWR']],
        [['LELB', 'LWR']],
        [['LELB', 'LWR']],
        ],
    'radius_l': [
        [['LELB', 'LWR']],
        [['LELB', 'LWR']],
        [['LELB', 'LWR']],
        ],
    }

#Set model path
pathModelFolder = os.getcwd()+'\\data\\Model\\'

#Set up logging.
logPath = os.path.join(pathModelFolder,'scaling.log')
if os.path.exists(logPath):
    os.remove(logPath)
    
#Remove all handlers associated with the root logger object.
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
logging.shutdown()
logging.basicConfig(filename = logPath,format='%(message)s',
                    level = logging.INFO)
osim.Logger.setLevelString('error')

#Load in the two models
laiUhlrich = osim.Model(pathModelFolder+'LaiUhlrich2022.osim')
laiUhlrichState = laiUhlrich.initSystem()
jaScaled = osim.Model(pathModelFolder+'JA1_SCALED_Osim40_Muscles.osim')
jaScaledState = jaScaled.initSystem()

#Calculate the scaling factors between marker pairs
scaleFactors = {}

#Loop through bodies
for body in scaleMarkerPairs.keys():
    
    #Add empty array in scale factors dictionary
    scaleFactors[body] = np.zeros((3,))
    
    #Loop through three axes
    for axInd in range(3):
        
        #Calculate distances for each model for the current body and axis
        #JA model
        jaModelDist = np.array([np.linalg.norm(jaScaled.updMarkerSet().get(scaleMarkerPairs[body][axInd][markerPairInd][0]).getLocationInGround(jaScaledState).to_numpy() - \
                                               jaScaled.updMarkerSet().get(scaleMarkerPairs[body][axInd][markerPairInd][1]).getLocationInGround(jaScaledState).to_numpy()) for \
                                markerPairInd in range(len(scaleMarkerPairs[body][axInd]))])
        #LaiUhlrich model
        laiModelDist = np.array([np.linalg.norm(laiUhlrich.updMarkerSet().get(scaleMarkerPairs[body][axInd][markerPairInd][0]).getLocationInGround(laiUhlrichState).to_numpy() - \
                                                laiUhlrich.updMarkerSet().get(scaleMarkerPairs[body][axInd][markerPairInd][1]).getLocationInGround(laiUhlrichState).to_numpy()) for \
                                 markerPairInd in range(len(scaleMarkerPairs[body][axInd]))])
            
        #Calculate scale factor average and append to scaling dictionary
        scaleFactors[body][axInd] = np.mean(jaModelDist / laiModelDist)
        
#Create scale tool
scaleTool = osim.ScaleTool()

#Set model mass
jaMass = np.sum([jaScaled.updBodySet().get(bodyInd).getMass() for bodyInd in range(jaScaled.updBodySet().getSize())])
scaleTool.setSubjectMass(jaMass)

#Set unscaled model to LaiUlrich
scaleTool.getGenericModelMaker().setModelFileName('data\Model\\LaiUhlrich2022.osim')

#Set scaling order to use manual measurements
scalingOrder = osim.ArrayStr()
scalingOrder.append('manualScale')
scaleTool.getModelScaler().setScalingOrder(scalingOrder)

#Set measurement scales in scale set
for body in scaleFactors.keys():
    
    #Create scale
    scale = osim.Scale()
    
    #Set parameters
    scale.setScaleFactors(osim.Vec3(scaleFactors[body]))
    scale.setSegmentName(body)
    
    #Append to scale set in tool
    scaleTool.getModelScaler().getScaleSet().cloneAndAppend(scale)

#Set output model file
scaleTool.getModelScaler().setOutputModelFileName('data\Model\\LaiUhlrich2022_JA1_SCALED.osim')

#Set marker placer to false
scaleTool.getMarkerPlacer().setApply(False)

#Run scale tool
scaleTool.run()

#Read in model to make some updates
scaledModel = osim.Model('data\Model\\LaiUhlrich2022_JA1_SCALED.osim')

#Delete markerset as this isn't needed
scaledModel.updMarkerSet().clearAndDestroy()

#Update name
scaledModel.setName('LaiUhlrich2022_JA1_SCALED')

#Re-print to file
scaledModel.printToXML('data\Model\\LaiUhlrich2022_JA1_SCALED.osim')
    
    
# %% Set data

#Adjust kinematic datafile to invert knee angle to match with new model
#Convert to joint angle headers .mot file too

#Read in data
coordinatesData = osim.TimeSeriesTable(os.path.join('data', 'Kinematics', trial_name+'.sto'))
    
#Get column headers
colHeaders = coordinatesData.getColumnLabels()

#Rename to remove the joint coordinate state details
newColHeaders = []
for col in colHeaders:
    if col.endswith('/value'):
        #Split the string to the 3rd component which is the joint angle
        newColHeaders.append(col.split('/')[3])
    else:
        newColHeaders.append(col)
        
#Set new column headers in coordinates data
coordinatesData.setColumnLabels(newColHeaders)

#Invert knee angle data
knee_angle_r_inverted = coordinatesData.getDependentColumn('knee_angle_r').to_numpy() * -1
knee_angle_l_inverted = coordinatesData.getDependentColumn('knee_angle_l').to_numpy() * -1

#Set in time series table by extracting vector and then row vector
for ii in range(len(coordinatesData.getIndependentColumn())):
    #Right knee angle at index
    coordinatesData.getDependentColumn('knee_angle_r').row(ii).setTo(knee_angle_r_inverted[ii])
    #Left knee angle at index
    coordinatesData.getDependentColumn('knee_angle_l').row(ii).setTo(knee_angle_l_inverted[ii])
    
#Write to .mot file
osim.STOFileAdapter().write(coordinatesData, os.path.join('data', 'Kinematics', trial_name+'.mot'))
    
# %% Adjust muscle wrapping

"""

The code in this section comes from the adjustMuscleWrapping function included 
in the utilsProcessing script included with the opencap-processing example.

"""

#Define a function to set poses and get moment arms
def getMomentArms(model, poses, muscleName, coordinateForMomentArm):
    state = model.initSystem()
    coords = model.getCoordinateSet()
    muscleSet = model.getMuscles()
    coordForMA = coords.get(coordinateForMomentArm)
    momentArms = []
    for i, pose in enumerate(poses):        
        for coordVal in pose:
            coords.get(coordVal[0]).setValue(state,np.deg2rad(coordVal[1]))
        momentArms.append(
            muscleSet.get(muscleName).computeMomentArm(state,coordForMA))
        
    return momentArms

#Start adjusting muscle wrapping surfaces
print('Adjust muscle wrapping surfaces.')
    
#Set up logging.
logPath = os.path.join(pathModelFolder,'modelAdjustment.log')
if os.path.exists(logPath):
    os.remove(logPath)
    
#Remove all handlers associated with the root logger object.
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
logging.shutdown()
logging.basicConfig(filename = logPath,format='%(message)s',
                    level = logging.INFO)

#Load models.
osim.Logger.setLevelString('error')
unscaledModel = osim.Model(pathModelFolder+'LaiUhlrich2022.osim')
scaledModel = osim.Model(pathModelFolder+'LaiUhlrich2022_JA1_SCALED.osim')
pathOutputModel = pathModelFolder+'LaiUhlrich2022_JA1_SCALED_adjusted.osim'
scaledBodySet = scaledModel.getBodySet()

#Poses that often cause problems.
pose_gmax = [
    [['hip_flexion_r',90],['hip_adduction_r',-26], ['hip_rotation_r',40]]
    ]
coord_gmax = 'hip_flexion_r'

#Generic model doesn't wrap beyond 32deg abd.
pose_hipFlexors = [
    [['hip_flexion_r',-30],['hip_adduction_r',-32],['hip_rotation_r',-36]],
    [['hip_flexion_r',-30],['hip_adduction_r',-50],['hip_rotation_r',0]],
    [['hip_flexion_r',-30],['hip_adduction_r',30],['hip_rotation_r',0]]] 
coord_hipFlexors = 'hip_flexion_r'

#Gmax1 - shrink wrap cyl radius.
momentArmsGmax_unscaled = getMomentArms(
    unscaledModel,pose_gmax,'glmax1_r',coord_gmax)
momentArmsGmax_scaled = getMomentArms(
    scaledModel,pose_gmax,'glmax1_r',coord_gmax)

#Get wrapping surface.
pelvis = scaledBodySet.get('pelvis')
gmaxWrap = osim.WrapCylinder.safeDownCast(
    pelvis.getWrapObjectSet().get('Gmax1_at_pelvis_r'))
radius = gmaxWrap.get_radius()
originalRadius = np.copy(radius)

#Iteratively adjust wrapping surface radius until desirable outcome is met
for iPose,(momentArmGmax_scaled,momentArmGmax_unscaled) in enumerate(zip(momentArmsGmax_scaled,momentArmsGmax_unscaled)): 
    if np.abs(momentArmGmax_scaled) < np.max([0.5* np.abs(momentArmGmax_unscaled), 0.008]): # This constant came from 100 scaled models
        originalBadMomentArm = np.copy(momentArmGmax_scaled)            
        while np.abs(momentArmGmax_scaled) <= np.abs(originalBadMomentArm) and radius > 0.002:
            gmaxWrap.set_radius(radius-0.002) 
            momentArmGmax_scaled = getMomentArms(scaledModel,pose_gmax,'glmax1_r',coord_gmax)[iPose]
            radius = gmaxWrap.get_radius()                
        if radius > 0.5*originalRadius:
            outputStr = '-For pose #{}, scaled gmax1 moment arm was {:.3f}. Unscaled was {:.3f}. Reduced R&L wrap radius from {:.3f} to {:.3f}, which increased the moment arm back to {:.3f}.'.format(
                          iPose, originalBadMomentArm,momentArmGmax_unscaled,
                          originalRadius,radius,momentArmGmax_scaled)
            print(outputStr)
            logging.info(outputStr)
            # Set the left side as well.
            osim.WrapCylinder.safeDownCast(pelvis.getWrapObjectSet().get('Gmax1_at_pelvis_l')).set_radius(radius)        
        else:
            outputStr = '-For pose #{}, couldn''t restore glmax1 moment arm by shrinking radius by 1/2. Model unchanged.'.format(iPose)
            print(outputStr)
            logging.info(outputStr)
            gmaxWrap.set_radius(float(originalRadius))        
        scaledModel.initSystem()       
    else:
        outputStr = '-For pose #{}, scaled gmax1 moment arm was {:.3f}. Unscaled was {:.3f}. No adjustments made.'.format(
                     iPose,np.abs(momentArmGmax_scaled),np.abs(momentArmGmax_unscaled))
        print(outputStr)
        logging.info(outputStr)

#Iliacus - change path points to engage wrap cylinder.
momentArms_unscaled = getMomentArms(
    unscaledModel,pose_hipFlexors,'iliacus_r',coord_hipFlexors)
momentArms_scaled = getMomentArms(
    scaledModel,pose_hipFlexors,'iliacus_r',coord_hipFlexors)    

#Get path point locations.
muscle = scaledModel.getMuscles().get('iliacus_r')
pathPoints = muscle.get_GeometryPath().getPathPointSet()
point1 = osim.PathPoint.safeDownCast(pathPoints.get(1))
loc1Vec = point1.get_location()
point2 = osim.PathPoint.safeDownCast(pathPoints.get(2))
loc2Vec = point2.get_location()    
original_loc1 = [loc1Vec[i] for i in range(3)]
original_loc2 = [loc2Vec[i] for i in range(3)]

#Get wrap cyl.        
wrapCyl = osim.WrapCylinder.safeDownCast(
    pelvis.getWrapObjectSet().get('IL_at_brim_r'))
radius = wrapCyl.get_radius()
originalRadius = np.copy(radius)
previousRadius = np.copy(radius)

#Iteratively adjust path points and wrapping surface radius until desirable outcome is met
for iPose,(momentArm_scaled,momentArm_unscaled) in enumerate(zip(momentArms_scaled,momentArms_unscaled)):
    if np.abs(momentArm_scaled) < np.max([0.7* np.abs(momentArm_unscaled) , 0.015]):             
        # Get path point locations.
        muscle = scaledModel.getMuscles().get('iliacus_r')
        pathPoints = muscle.get_GeometryPath().getPathPointSet()
        point1 = osim.PathPoint.safeDownCast(pathPoints.get(1))
        loc1Vec = point1.get_location()
        point2 = osim.PathPoint.safeDownCast(pathPoints.get(2))
        loc2Vec = point2.get_location()            
        originalBadMomentArm = np.copy(momentArm_scaled)                           
        while np.abs(momentArm_scaled) <= np.max([0.7* np.abs(momentArm_unscaled) , 0.015]) and (np.abs(loc1Vec[0]-original_loc1[0]) < 0.015 and np.abs(loc2Vec[1]-original_loc2[1]) <0.015):
            loc1Vec[0] += 0.002 # Move the 1st (pelvis) path point forward
            loc2Vec[1] -= 0.002 # move the 2nd (femur) path point down
            point1.set_location(loc1Vec)
            point2.set_location(loc2Vec)        
            momentArm_scaled = getMomentArms(scaledModel,pose_hipFlexors,'iliacus_r',coord_hipFlexors)[iPose]          
        while np.abs(momentArm_scaled) <= np.max([0.7* np.abs(momentArm_unscaled) , 0.015]) and radius>0.7*originalRadius: # above approach did not succeed, drop the cyl radius some
            wrapCyl.set_radius(radius-0.002) 
            momentArm_scaled = getMomentArms(scaledModel,pose_hipFlexors,'iliacus_r',coord_hipFlexors)[iPose]
            pelvis = scaledBodySet.get('pelvis')
            radius = wrapCyl.get_radius()
        if np.abs(momentArm_scaled) > np.max([0.7* np.abs(momentArm_unscaled) , 0.015]): # succeeded
            # Set the left side as well.
            muscle = scaledModel.getMuscles().get('iliacus_l')
            pathPoints = muscle.get_GeometryPath().getPathPointSet()        
            point1 = osim.PathPoint.safeDownCast(pathPoints.get(1))
            loc1Vec_l = point1.get_location()
            loc1Vec_l[0] = loc1Vec[0]
            point1.set_location(loc1Vec_l)                
            point2 = osim.PathPoint.safeDownCast(pathPoints.get(2))
            loc2Vec_l = point2.get_location()
            loc2Vec_l[1] = loc2Vec[1]
            point2.set_location(loc2Vec_l)                
            if radius<previousRadius:
                radiusStr = ', and after moving points by 1.5±0.2cm wasn''t enough, reduced R&L iliacus wrap radius from {:.3f} to {:.3f}'.format(
                originalRadius,radius)
                # set the left side as well.
                osim.WrapCylinder.safeDownCast(pelvis.getWrapObjectSet().get('IL_at_brim_l')).set_radius(radius)
            else:
                radiusStr = ''
            previousRadius = np.copy(radius)    
            outputStr = '-For pose #{}, moved iliacus pelvis path point xpos forward from {:.3f} to {:.3f}, and femur iliacus path point ypos down from {:.3f} to {:.3f}'.format(
                iPose,original_loc1[0],loc1Vec[0],original_loc2[1],loc2Vec[1]) + radiusStr + '. Restored moment arm from {:.3f} to {:.3f}.'.format(
                  originalBadMomentArm,momentArm_scaled)
            print(outputStr)
            logging.info(outputStr)
        else:
            outputStr = '-For pose #{}, couldn''t restore iliacus moment arm by moving path points by 2cm. Model unchanged.'.format(iPose)
            print(outputStr)
            logging.info(outputStr)                
            point1.set_location(original_loc1)
            point2.set_location(original_loc2)            
        scaledModel.initSystem()
    else:
        outputStr = '-For pose #{}, scaled iliacus moment arm was {:.3f}. Unscaled was {:.3f}. No adjustments made.'.format(
              iPose,np.abs(momentArm_scaled),np.abs(momentArm_unscaled))
        print(outputStr)
        logging.info(outputStr)

#Psoas - change path points to engage wrap cylinder.
momentArms_unscaled = getMomentArms(
    unscaledModel,pose_hipFlexors,'psoas_r',coord_hipFlexors)
momentArms_scaled = getMomentArms(
    scaledModel,pose_hipFlexors,'psoas_r',coord_hipFlexors)   
 
#Get path point locations 
muscle = scaledModel.getMuscles().get('psoas_r')
pathPoints = muscle.get_GeometryPath().getPathPointSet()
point1 = osim.PathPoint.safeDownCast(pathPoints.get(1))
loc1Vec = point1.get_location()
point2 = osim.PathPoint.safeDownCast(pathPoints.get(2))
loc2Vec = point2.get_location()    
original_loc1 = [loc1Vec[i] for i in range(3)]
original_loc2 = [loc2Vec[i] for i in range(3)]

#Get wrap cyl         
wrapCyl = osim.WrapCylinder.safeDownCast(
    pelvis.getWrapObjectSet().get('PS_at_brim_r'))
radius = wrapCyl.get_radius()
originalRadius = np.copy(radius)
previousRadius = np.copy(radius)

#Iteratively adjust path points and wrapping surface radius until desirable outcome is met
for iPose,(momentArm_scaled,momentArm_unscaled) in enumerate(zip(momentArms_scaled,momentArms_unscaled)):
    if np.abs(momentArm_scaled) < np.max([0.7* np.abs(momentArm_unscaled), 0.015]):            
        # Get path point locations.
        muscle = scaledModel.getMuscles().get('psoas_r')
        pathPoints = muscle.get_GeometryPath().getPathPointSet()
        point1 = osim.PathPoint.safeDownCast(pathPoints.get(1))
        loc1Vec = point1.get_location()
        point2 = osim.PathPoint.safeDownCast(pathPoints.get(2))
        loc2Vec = point2.get_location()
        originalBadMomentArm = np.copy(momentArm_scaled)               
        while np.abs(momentArm_scaled) <= np.max([0.7* np.abs(momentArm_unscaled), 0.015]) and (np.abs(loc1Vec[0]-original_loc1[0]) < 0.015 and np.abs(loc2Vec[1]-original_loc2[1]) < 0.015):
            loc1Vec[0] += 0.002 # Move the 1st (pelvis) path point forward
            loc2Vec[1] -= 0.002 # move the 2nd (femur) path point down
            point1.set_location(loc1Vec)
            point2.set_location(loc2Vec)        
            momentArm_scaled = getMomentArms(scaledModel,pose_hipFlexors,'psoas_r',coord_hipFlexors)[iPose]            
        while np.abs(momentArm_scaled) <= np.max([0.7* np.abs(momentArm_unscaled) , 0.015]) and radius>0.7*originalRadius: #above approach did not succeed, drop the cyl radius some
            wrapCyl.set_radius(radius-0.002) 
            momentArm_scaled = getMomentArms(scaledModel,pose_hipFlexors,'psoas_r',coord_hipFlexors)[iPose]
            pelvis = scaledBodySet.get('pelvis')
            radius = wrapCyl.get_radius()
        if np.abs(momentArm_scaled) > np.max([0.7* np.abs(momentArm_unscaled) , 0.015]): # succeeded
            # set the left side as well.
            muscle = scaledModel.getMuscles().get('psoas_l')
            pathPoints = muscle.get_GeometryPath().getPathPointSet()        
            point1 = osim.PathPoint.safeDownCast(pathPoints.get(1))
            loc1Vec_l = point1.get_location()
            loc1Vec_l[0] = loc1Vec[0]
            point1.set_location(loc1Vec_l)
            point2 = osim.PathPoint.safeDownCast(pathPoints.get(2))
            loc2Vec_l = point2.get_location()
            loc2Vec_l[1] = loc2Vec[1]
            point2.set_location(loc2Vec_l)
            if radius<previousRadius:
                radiusStr = ', and after moving points by 1.5±0.2cm wasn''t enough, reduced R&L psoas wrap radius from {:.3f} to {:.3f}'.format(
                originalRadius,radius)
                # set the left side as well.
                osim.WrapCylinder.safeDownCast(pelvis.getWrapObjectSet().get('PS_at_brim_l')).set_radius(radius)
            else:
                radiusStr = ''
            previousRadius = np.copy(radius)   
            outputStr = '-For pose #{}, moved psoas pelvis path point xpos forward from {:.3f} to {:.3f}, and femur psoas path point ypos down from {:.3f} to {:.3f}'.format(
                iPose,original_loc1[0],loc1Vec[0],original_loc2[1],loc2Vec[1]) + radiusStr + '. Restored moment arm from {:.3f} to {:.3f}.'.format(
                  originalBadMomentArm,momentArm_scaled)
            print(outputStr)
            logging.info(outputStr)                   
        else:
            outputStr = '-For pose #{}, couldn''t restore psoas moment arm by moving path points by 2cm. Model unchanged.'.format(iPose)
            print(outputStr)
            logging.info(outputStr)                
            point1.set_location(osim.Vec3(original_loc1))
            point2.set_location(osim.Vec3(original_loc2))          
        scaledModel.initSystem()           
    else:
        outputStr = '-For pose #{}, scaled psoas moment arm was {:.3f}. Unscaled was {:.3f}. No adjustements made.'.format(
              iPose,np.abs(momentArm_scaled),np.abs(momentArm_unscaled))
        print(outputStr)
        logging.info(outputStr)

#Reprint model to file
scaledModel.printToXML(pathOutputModel)
logging.shutdown()

# %% Add foot ground contacts

"""

The code in this section comes from the generateModelWithContacts function
included in the utilsProcessing script included with the opencap-processing
example.

"""

#Set output model path
adjustedModelFile = pathModelFolder+'LaiUhlrich2022_JA1_SCALED_adjusted.osim'
pathOutputModel = pathModelFolder+'LaiUhlrich2022_JA1_SCALED_adjusted_contacts.osim'

#Settings for patella mass
setPatellaMassToZero = True

#Start adding foot-ground contacts
print('Add foot-ground contacts.')
        
#Add contact spheres to the scaled model.
#The parameters of the foot-ground contacts are based on previous work. We
#scale the contact sphere locations based on foot dimensions.

#Set the reference contact spheres
#Note that these come from the original example script
#TODO: do these need to be moved?
reference_contact_spheres = {
    "s1_r": {"radius": 0.032, "location": np.array([0.0019011578840796601,   -0.01,  -0.00382630379623308]), "orientation": np.array([0, 0, 0]), "socket_frame": "calcn_r"},
    "s2_r": {"radius": 0.032, "location": np.array([0.14838639994206301,     -0.01,  -0.028713422052654002]), "orientation": np.array([0, 0, 0]), "socket_frame": "calcn_r"},
    "s3_r": {"radius": 0.032, "location": np.array([0.13300117060705099,     -0.01,  0.051636247344956601]), "orientation": np.array([0, 0, 0]), "socket_frame": "calcn_r"},
    "s4_r": {"radius": 0.032, "location": np.array([0.066234666199163503,    -0.01,  0.026364160674169801]), "orientation": np.array([0, 0, 0]), "socket_frame": "calcn_r"},
    "s5_r": {"radius": 0.032, "location": np.array([0.059999999999999998,    -0.01,  -0.018760308461917698]), "orientation": np.array([0, 0, 0]), "socket_frame": "toes_r" },
    "s6_r": {"radius": 0.032, "location": np.array([0.044999999999999998,    -0.01,  0.061856956754965199]), "orientation": np.array([0, 0, 0]), "socket_frame": "toes_r" },
    "s1_l": {"radius": 0.032, "location": np.array([0.0019011578840796601,   -0.01,  0.00382630379623308]), "orientation": np.array([0, 0, 0]), "socket_frame": "calcn_l"},
    "s2_l": {"radius": 0.032, "location": np.array([0.14838639994206301,     -0.01,  0.028713422052654002]), "orientation": np.array([0, 0, 0]), "socket_frame": "calcn_l"},
    "s3_l": {"radius": 0.032, "location": np.array([0.13300117060705099,     -0.01,  -0.051636247344956601]), "orientation": np.array([0, 0, 0]), "socket_frame": "calcn_l"},
    "s4_l": {"radius": 0.032, "location": np.array([0.066234666199163503,    -0.01,  -0.026364160674169801]), "orientation": np.array([0, 0, 0]), "socket_frame": "calcn_l"},
    "s5_l": {"radius": 0.032, "location": np.array([0.059999999999999998,    -0.01,  0.018760308461917698]), "orientation": np.array([0, 0, 0]), "socket_frame": "toes_l" },
    "s6_l": {"radius": 0.032, "location": np.array([0.044999999999999998,    -0.01,  -0.061856956754965199]), "orientation": np.array([0, 0, 0]), "socket_frame": "toes_l" }}      

#Set the reference scale factors
#Note that these come from the original example script
#TODO: do these need to be changed?
reference_scale_factors = {"calcn_r": np.array([0.91392399999999996, 0.91392399999999996, 0.91392399999999996]),
                           "toes_r":  np.array([0.91392399999999996, 0.91392399999999996, 0.91392399999999996]),
                           "calcn_l": np.array([0.91392399999999996, 0.91392399999999996, 0.91392399999999996]),
                           "toes_l":  np.array([0.91392399999999996, 0.91392399999999996, 0.91392399999999996])}

#Set the contact half space
reference_contact_half_space = {"name": "floor", "location": np.array([0, 0, 0]),"orientation": np.array([0, 0, -np.pi/2]), "frame": "ground"}

#Set contact sphere parameters
stiffness = 1000000
dissipation = 2.0
static_friction = 0.8
dynamic_friction = 0.8
viscous_friction = 0.5
transition_velocity = 0.2

#Add contact spheres and SmoothSphereHalfSpaceForces.
osim.Logger.setLevelString('error')
model = osim.Model(adjustedModelFile)
bodySet = model.get_BodySet()

#ContactHalfSpace.
if reference_contact_half_space["frame"] == "ground":
    contact_half_space_frame = model.get_ground()
else:
    raise ValueError('Not yet supported.')    
contactHalfSpace = osim.ContactHalfSpace(
    osim.Vec3(reference_contact_half_space["location"]),
    osim.Vec3(reference_contact_half_space["orientation"]),
    contact_half_space_frame, reference_contact_half_space["name"])
contactHalfSpace.connectSocket_frame(contact_half_space_frame)
model.addContactGeometry(contactHalfSpace)

#ContactSpheres and SmoothSphereHalfSpaceForces.
for ref_contact_sphere in reference_contact_spheres:    
    # ContactSpheres.
    body = bodySet.get(reference_contact_spheres[ref_contact_sphere]["socket_frame"])
    # Scale location based on attached_geometry scale_factors.      
    # We don't scale the y_position.
    attached_geometry = body.get_attached_geometry(0)
    c_scale_factors = attached_geometry.get_scale_factors().to_numpy() 
    c_ref_scale_factors = reference_scale_factors[reference_contact_spheres[ref_contact_sphere]["socket_frame"]]
    scale_factors = c_ref_scale_factors / c_scale_factors        
    scale_factors[1] = 1        
    scaled_location = reference_contact_spheres[ref_contact_sphere]["location"] / scale_factors
    c_contactSphere = osim.ContactSphere(
        reference_contact_spheres[ref_contact_sphere]["radius"],
        osim.Vec3(scaled_location), body, ref_contact_sphere)
    c_contactSphere.connectSocket_frame(body)
    model.addContactGeometry(c_contactSphere)
    
    # SmoothSphereHalfSpaceForces.
    SmoothSphereHalfSpaceForce = osim.SmoothSphereHalfSpaceForce(
        "SmoothSphereHalfSpaceForce_" + ref_contact_sphere, 
        c_contactSphere, contactHalfSpace)
    SmoothSphereHalfSpaceForce.set_stiffness(stiffness)
    SmoothSphereHalfSpaceForce.set_dissipation(dissipation)
    SmoothSphereHalfSpaceForce.set_static_friction(static_friction)
    SmoothSphereHalfSpaceForce.set_dynamic_friction(dynamic_friction)
    SmoothSphereHalfSpaceForce.set_viscous_friction(viscous_friction)
    SmoothSphereHalfSpaceForce.set_transition_velocity(transition_velocity)        
    SmoothSphereHalfSpaceForce.connectSocket_half_space(contactHalfSpace)
    SmoothSphereHalfSpaceForce.connectSocket_sphere(c_contactSphere)
    model.addForce(SmoothSphereHalfSpaceForce)

# We do not use the patella in the dynamic simulations. The reason is that
# the patella only matters for the muscle-tendon lengths and moment arms,
# but since we approximate those with polynomials, the patella is useless.
# We therefore remove it, since otherwise we would have to deal with
# kinematic constraints that would make things unecessarily complicated.
# We remove it when building the external function, and here we set its
# mass to zero such that we can make an apple-to-apple comparison when
# checking that the outputs from the external function match the results
# from ID ran with the model (with a mass set to 0, the patella will not
# influence ID).
if setPatellaMassToZero:
    for i in range(bodySet.getSize()):        
        c_body = bodySet.get(i)
        c_body_name = c_body.getName()            
        if (c_body_name == 'patella_l' or c_body_name == 'patella_r'):
            c_body.set_mass(0.)
            c_body.set_inertia(osim.Vec6(0))

#Finalize model and print to file
model.finalizeConnections
model.initSystem()
model.printToXML(pathOutputModel)

# %% Generate external functions

"""

The code in this section comes from the generateExternalFunction function
included in the utilsOpenSimAD script included with the opencap-processing
example. This is a huge function, and hence would make sense to include in an
external script. Most of this comes from building the external function elements
into the file.

"""

#Set whether using a treadmill
treadmill = False

#Set overwrite to false to avoid re-doing
overwrite = False

#Set model file
outputModelFileName = pathModelFolder+'LaiUhlrich2022_JA1_SCALED_adjusted_contacts.osim'

#Set external functions folder
pathOutputExternalFunctionFolder = os.path.join(pathModelFolder, 'ExternalFunction')
os.makedirs(pathOutputExternalFunctionFolder, exist_ok = True)

#Set external function name
externalFunctionName = 'F'
#Add treadmill label if needed
if treadmill:
    externalFunctionName += '_treadmill'    

#Set outut file names
pathOutputFile = os.path.join(pathOutputExternalFunctionFolder, externalFunctionName + '.cpp')
pathOutputMap = os.path.join(pathOutputExternalFunctionFolder, externalFunctionName + '_map.npy')

#Set external function file extension
if platform.system() == 'Windows':
    ext_F = '.dll'
elif platform.system() == 'Darwin':
    ext_F = '.dylib'
elif platform.system() == 'Linux':
    ext_F = '.so'
else:
    raise ValueError("Platform not supported.")
    
#Set the full library output path
pathOutputDll = os.path.join(pathOutputExternalFunctionFolder, externalFunctionName + ext_F)

# if (overwrite is False and os.path.exists(pathOutputFile) and 
#     os.path.exists(pathOutputMap) and os.path.exists(pathOutputDll)):
#     return      
# else:

#Get ready to generate external function
print('Generate external function to leverage automatic differentiation.')

#Generate external Function (.cpp file)

#Set-up the logger
osim.Logger.setLevelString('error')

#Read in model and details
model = osim.Model(outputModelFileName)
model.initSystem()
bodySet = model.getBodySet()
jointSet = model.get_JointSet()
nJoints = jointSet.getSize()
geometrySet = model.get_ContactGeometrySet()
forceSet = model.get_ForceSet()
coordinateSet = model.getCoordinateSet()
nCoordinates = coordinateSet.getSize()

#Create a coordinates list
coordinates = []
for coor in range(nCoordinates):
    coordinates.append(coordinateSet.get(coor).getName())

#Remove patellofemoral joint
sides = ['r', 'l']
for side in sides:
    # We do not include the coordinates from the patellofemoral joints,
    # since they only influence muscle paths, which we approximate using
    # polynomials.
    if 'knee_angle_{}_beta'.format(side) in coordinates:
        nCoordinates -= 1
        nJoints -= 1

#Calculate number of bodies minus patella
nBodies = 0
for i in range(bodySet.getSize()):        
    c_body = bodySet.get(i)
    c_body_name = c_body.getName()  
    if (c_body_name == 'patella_l' or c_body_name == 'patella_r'):
        continue
    nBodies += 1

#Calculate number of contacts
nContacts = 0
for i in range(forceSet.getSize()):        
    c_force_elt = forceSet.get(i)        
    if c_force_elt.getConcreteClassName() == "SmoothSphereHalfSpaceForce":  
        nContacts += 1

#Write external function file
with open(pathOutputFile, "w") as f:    
    
    #OpenSim libraries    
    # TODO: only include those that are necessary (model-specific).
    f.write('#include <OpenSim/Simulation/Model/Model.h>\n')
    f.write('#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>\n')
    f.write('#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>\n')
    f.write('#include <OpenSim/Simulation/SimbodyEngine/Joint.h>\n')
    f.write('#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>\n')
    f.write('#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>\n')
    if treadmill:
        f.write('#include <OpenSim/Simulation/SimbodyEngine/SliderJoint.h>\n')    
    f.write('#include <OpenSim/Common/LinearFunction.h>\n')
    f.write('#include <OpenSim/Common/PolynomialFunction.h>\n')
    f.write('#include <OpenSim/Common/MultiplierFunction.h>\n')
    f.write('#include <OpenSim/Common/Constant.h>\n')
    f.write('#include <OpenSim/Simulation/Model/SmoothSphereHalfSpaceForce.h>\n')
    f.write('#include <OpenSim/Simulation/SimulationUtilities.h>\n')
    f.write('#include "SimTKcommon/internal/recorder.h"\n\n')
    
    f.write('#include <iostream>\n')
    f.write('#include <iterator>\n')
    f.write('#include <random>\n')
    f.write('#include <cassert>\n')
    f.write('#include <algorithm>\n')
    f.write('#include <vector>\n')
    f.write('#include <fstream>\n\n')
    
    f.write('using namespace SimTK;\n')
    f.write('using namespace OpenSim;\n\n')

    if treadmill:
        f.write('constexpr int n_in = 3; \n')
    else:
        f.write('constexpr int n_in = 2; \n')
    f.write('constexpr int n_out = 1; \n')
    
    #Coordinates
    f.write('constexpr int nCoordinates = %i; \n' % nCoordinates)
    f.write('constexpr int NX = nCoordinates*2; \n')
    f.write('constexpr int NU = nCoordinates; \n\n')
    if treadmill:
        nCoordinates_treadmill = nCoordinates + 1
        f.write('constexpr int nCoordinates_treadmill = %i; \n' % nCoordinates_treadmill)
        f.write('constexpr int NX_treadmill = nCoordinates_treadmill*2; \n')
        f.write('constexpr int NU_treadmill = nCoordinates_treadmill; \n\n')

    f.write('template<typename T> \n')
    f.write('T value(const Recorder& e) { return e; }; \n')
    f.write('template<> \n')
    f.write('double value(const Recorder& e) { return e.getValue(); }; \n\n')
    
    f.write('template<typename T>\n')
    f.write('int F_generic(const T** arg, T** res) {\n\n')
    
    #Model
    f.write('\t// Definition of model.\n')
    f.write('\tOpenSim::Model* model;\n')
    f.write('\tmodel = new OpenSim::Model();\n\n')
    
    #Bodies
    f.write('\t// Definition of bodies.\n')
    for i in range(bodySet.getSize()):        
        c_body = bodySet.get(i)
        c_body_name = c_body.getName()            
        if (c_body_name == 'patella_l' or c_body_name == 'patella_r'):
            continue            
        c_body_mass = c_body.get_mass()
        c_body_mass_center = c_body.get_mass_center().to_numpy()
        c_body_inertia = c_body.get_inertia()
        c_body_inertia_vec3 = np.array([c_body_inertia.get(0), c_body_inertia.get(1), c_body_inertia.get(2)])        
        f.write('\tOpenSim::Body* %s;\n' % c_body_name)
        f.write('\t%s = new OpenSim::Body(\"%s\", %.20f, Vec3(%.20f, %.20f, %.20f), Inertia(%.20f, %.20f, %.20f, 0., 0., 0.));\n' % (c_body_name, c_body_name, c_body_mass, c_body_mass_center[0], c_body_mass_center[1], c_body_mass_center[2], c_body_inertia_vec3[0], c_body_inertia_vec3[1], c_body_inertia_vec3[2]))
        f.write('\tmodel->addBody(%s);\n' % (c_body_name))
        f.write('\n')
    if treadmill:
        f.write('\tOpenSim::Body* treadmill;\n')
        f.write('\ttreadmill = new OpenSim::Body("treadmill", 1., Vec3(0), Inertia(1,1,1,0,0,0));\n')
        f.write('\tmodel->addBody(treadmill);\n')
        f.write('\n')
    
    #Joints
    f.write('\t// Definition of joints.\n')
    for i in range(jointSet.getSize()): 
        c_joint = jointSet.get(i)
        c_joint_type = c_joint.getConcreteClassName()
        
        c_joint_name = c_joint.getName()
        if (c_joint_name == 'patellofemoral_l' or 
            c_joint_name == 'patellofemoral_r'):
            continue
        
        parent_frame = c_joint.get_frames(0)
        parent_frame_name = parent_frame.getParentFrame().getName()
        parent_frame_trans = parent_frame.get_translation().to_numpy()
        parent_frame_or = parent_frame.get_orientation().to_numpy()
        
        child_frame = c_joint.get_frames(1)
        child_frame_name = child_frame.getParentFrame().getName()
        child_frame_trans = child_frame.get_translation().to_numpy()
        child_frame_or = child_frame.get_orientation().to_numpy()
        
        #Custom joints
        if c_joint_type == "CustomJoint":
            
            f.write('\tSpatialTransform st_%s;\n' % c_joint.getName())                
            cObj = osim.CustomJoint.safeDownCast(c_joint)    
            spatialtransform = cObj.get_SpatialTransform()
            
            # Transform axis.
            # Rotation 1
            rot1 = spatialtransform.get_rotation1()
            rot1_axis = rot1.get_axis().to_numpy()
            rot1_f = rot1.get_function()
            coord = 0
            if rot1_f.getConcreteClassName() == 'LinearFunction':  
                rot1_f_obj = osim.LinearFunction.safeDownCast(rot1_f)                          
                rot1_f_slope = rot1_f_obj.getSlope()
                rot1_f_intercept = rot1_f_obj.getIntercept()                
                c_coord = c_joint.get_coordinates(coord)
                c_coord_name = c_coord.getName()
                f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                f.write('\tst_%s[%i].setFunction(new LinearFunction(%.4f, %.4f));\n' % (c_joint.getName(), coord, rot1_f_slope, rot1_f_intercept))                
            elif rot1_f.getConcreteClassName() == 'PolynomialFunction':
                f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                rot1_f_obj = osim.PolynomialFunction.safeDownCast(rot1_f)                
                rot1_f_coeffs = rot1_f_obj.getCoefficients().to_numpy()
                c_nCoeffs = rot1_f_coeffs.shape[0]                
                if c_nCoeffs == 2:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot1_f_coeffs[0], rot1_f_coeffs[1]))
                elif c_nCoeffs == 3:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot1_f_coeffs[0], rot1_f_coeffs[1], rot1_f_coeffs[2]))
                elif c_nCoeffs == 4:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot1_f_coeffs[0], rot1_f_coeffs[1], rot1_f_coeffs[2], rot1_f_coeffs[3]))  
                elif c_nCoeffs == 5:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot1_f_coeffs[0], rot1_f_coeffs[1], rot1_f_coeffs[2], rot1_f_coeffs[3], rot1_f_coeffs[4]))                    
                else:
                    raise ValueError("TODO")
                f.write('\tVector st_%s_%i_coeffs_vec(%i); \n' % (c_joint.getName(), coord, c_nCoeffs))
                f.write('\tfor (int i = 0; i < %i; ++i) st_%s_%i_coeffs_vec[i] = st_%s_%i_coeffs[i]; \n' % (c_nCoeffs, c_joint.getName(), coord, c_joint.getName(), coord))
                f.write('\tst_%s[%i].setFunction(new PolynomialFunction(st_%s_%i_coeffs_vec));\n' % (c_joint.getName(), coord, c_joint.getName(), coord))
            elif rot1_f.getConcreteClassName() == 'MultiplierFunction':
                rot1_f_obj = osim.MultiplierFunction.safeDownCast(rot1_f)
                rot1_f_obj_scale = rot1_f_obj.getScale()
                rot1_f_obj_f = rot1_f_obj.getFunction()
                rot1_f_obj_f_name = rot1_f_obj_f.getConcreteClassName()
                if rot1_f_obj_f_name == 'Constant':
                    rot1_f_obj_f_obj = osim.Constant.safeDownCast(rot1_f_obj_f)
                    rot1_f_obj_f_obj_value = rot1_f_obj_f_obj.getValue()
                    f.write('\tst_%s[%i].setFunction(new MultiplierFunction(new Constant(%.20f), %.20f));\n' % (c_joint.getName(), coord, rot1_f_obj_f_obj_value, rot1_f_obj_scale))
                elif rot1_f_obj_f_name == 'PolynomialFunction':
                    f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                    rot1_f_obj_f_obj = osim.PolynomialFunction.safeDownCast(rot1_f_obj_f)
                    rot1_f_obj_f_coeffs = rot1_f_obj_f_obj.getCoefficients().to_numpy()
                    c_nCoeffs = rot1_f_obj_f_coeffs.shape[0]
                    if c_nCoeffs == 2:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot1_f_obj_f_coeffs[0], rot1_f_obj_f_coeffs[1]))
                    elif c_nCoeffs == 3:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot1_f_obj_f_coeffs[0], rot1_f_obj_f_coeffs[1], rot1_f_obj_f_coeffs[2]))
                    elif c_nCoeffs == 4:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot1_f_obj_f_coeffs[0], rot1_f_obj_f_coeffs[1], rot1_f_obj_f_coeffs[2], rot1_f_obj_f_coeffs[3]))  
                    elif c_nCoeffs == 5:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot1_f_obj_f_coeffs[0], rot1_f_obj_f_coeffs[1], rot1_f_obj_f_coeffs[2], rot1_f_obj_f_coeffs[3], rot1_f_obj_f_coeffs[4]))                    
                    else:
                        raise ValueError("TODO")
                    f.write('\tVector st_%s_%i_coeffs_vec(%i); \n' % (c_joint.getName(), coord, c_nCoeffs))
                    f.write('\tfor (int i = 0; i < %i; ++i) st_%s_%i_coeffs_vec[i] = st_%s_%i_coeffs[i]; \n' % (c_nCoeffs, c_joint.getName(), coord, c_joint.getName(), coord))
                    f.write('\tst_%s[%i].setFunction(new MultiplierFunction(new PolynomialFunction(st_%s_%i_coeffs_vec), %.20f));\n' % (c_joint.getName(), coord, c_joint.getName(), coord, rot1_f_obj_scale))
                else:
                    raise ValueError("Not supported")
            elif rot1_f.getConcreteClassName() == 'Constant':
                rot1_f_obj = osim.Constant.safeDownCast(rot1_f)
                rot1_f_obj_value = rot1_f_obj.getValue()
                f.write('\tst_%s[%i].setFunction(new Constant(%.20f));\n' % (c_joint.getName(), coord, rot1_f_obj_value))
            else:
                raise ValueError("Not supported")
            f.write('\tst_%s[%i].setAxis(Vec3(%.20f, %.20f, %.20f));\n' % (c_joint.getName(), coord, rot1_axis[0], rot1_axis[1], rot1_axis[2]))
            
            #Rotation 2
            rot2 = spatialtransform.get_rotation2()
            rot2_axis = rot2.get_axis().to_numpy()
            rot2_f = rot2.get_function()
            coord = 1
            if rot2_f.getConcreteClassName() == 'LinearFunction':
                rot2_f_obj = osim.LinearFunction.safeDownCast(rot2_f)
                rot2_f_slope = rot2_f_obj.getSlope()
                rot2_f_intercept = rot2_f_obj.getIntercept()                
                c_coord = c_joint.get_coordinates(coord)
                c_coord_name = c_coord.getName()
                f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                f.write('\tst_%s[%i].setFunction(new LinearFunction(%.4f, %.4f));\n' % (c_joint.getName(), coord, rot2_f_slope, rot2_f_intercept))
            elif rot2_f.getConcreteClassName() == 'PolynomialFunction':
                f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                rot2_f_obj = osim.PolynomialFunction.safeDownCast(rot2_f)                
                rot2_f_coeffs = rot2_f_obj.getCoefficients().to_numpy()
                c_nCoeffs = rot2_f_coeffs.shape[0]                
                if c_nCoeffs == 2:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot2_f_coeffs[0], rot2_f_coeffs[1]))
                elif c_nCoeffs == 3:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot2_f_coeffs[0], rot2_f_coeffs[1], rot2_f_coeffs[2]))
                elif c_nCoeffs == 4:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot2_f_coeffs[0], rot2_f_coeffs[1], rot2_f_coeffs[2], rot2_f_coeffs[3]))  
                elif c_nCoeffs == 5:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot2_f_coeffs[0], rot2_f_coeffs[1], rot2_f_coeffs[2], rot2_f_coeffs[3], rot2_f_coeffs[4]))                    
                else:
                    raise ValueError("TODO")
                f.write('\tVector st_%s_%i_coeffs_vec(%i); \n' % (c_joint.getName(), coord, c_nCoeffs))
                f.write('\tfor (int i = 0; i < %i; ++i) st_%s_%i_coeffs_vec[i] = st_%s_%i_coeffs[i]; \n' % (c_nCoeffs, c_joint.getName(), coord, c_joint.getName(), coord))
                f.write('\tst_%s[%i].setFunction(new PolynomialFunction(st_%s_%i_coeffs_vec));\n' % (c_joint.getName(), coord, c_joint.getName(), coord))
            elif rot2_f.getConcreteClassName() == 'MultiplierFunction':
                rot2_f_obj = osim.MultiplierFunction.safeDownCast(rot2_f)
                rot2_f_obj_scale = rot2_f_obj.getScale()
                rot2_f_obj_f = rot2_f_obj.getFunction()
                rot2_f_obj_f_name = rot2_f_obj_f.getConcreteClassName()
                if rot2_f_obj_f_name == 'Constant':
                    rot2_f_obj_f_obj = osim.Constant.safeDownCast(rot2_f_obj_f)
                    rot2_f_obj_f_obj_value = rot2_f_obj_f_obj.getValue()
                    f.write('\tst_%s[%i].setFunction(new MultiplierFunction(new Constant(%.20f), %.20f));\n' % (c_joint.getName(), coord, rot2_f_obj_f_obj_value, rot2_f_obj_scale)) 
                elif rot2_f_obj_f_name == 'PolynomialFunction':
                    f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                    rot2_f_obj_f_obj = osim.PolynomialFunction.safeDownCast(rot2_f_obj_f)
                    rot2_f_obj_f_coeffs = rot2_f_obj_f_obj.getCoefficients().to_numpy()
                    c_nCoeffs = rot2_f_obj_f_coeffs.shape[0]
                    if c_nCoeffs == 2:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot2_f_obj_f_coeffs[0], rot2_f_obj_f_coeffs[1]))
                    elif c_nCoeffs == 3:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot2_f_obj_f_coeffs[0], rot2_f_obj_f_coeffs[1], rot2_f_obj_f_coeffs[2]))
                    elif c_nCoeffs == 4:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot2_f_obj_f_coeffs[0], rot2_f_obj_f_coeffs[1], rot2_f_obj_f_coeffs[2], rot2_f_obj_f_coeffs[3]))  
                    elif c_nCoeffs == 5:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot2_f_obj_f_coeffs[0], rot2_f_obj_f_coeffs[1], rot2_f_obj_f_coeffs[2], rot2_f_obj_f_coeffs[3], rot2_f_obj_f_coeffs[4]))                    
                    else:
                        raise ValueError("TODO")
                    f.write('\tVector st_%s_%i_coeffs_vec(%i); \n' % (c_joint.getName(), coord, c_nCoeffs))
                    f.write('\tfor (int i = 0; i < %i; ++i) st_%s_%i_coeffs_vec[i] = st_%s_%i_coeffs[i]; \n' % (c_nCoeffs, c_joint.getName(), coord, c_joint.getName(), coord))
                    f.write('\tst_%s[%i].setFunction(new MultiplierFunction(new PolynomialFunction(st_%s_%i_coeffs_vec), %.20f));\n' % (c_joint.getName(), coord, c_joint.getName(), coord, rot2_f_obj_scale))
                else:
                    raise ValueError("Not supported")
            elif rot2_f.getConcreteClassName() == 'Constant':
                rot2_f_obj = osim.Constant.safeDownCast(rot2_f)
                rot2_f_obj_value = rot2_f_obj.getValue()
                f.write('\tst_%s[%i].setFunction(new Constant(%.20f));\n' % (c_joint.getName(), coord, rot2_f_obj_value))
            else:
                raise ValueError("Not supported")
            f.write('\tst_%s[%i].setAxis(Vec3(%.20f, %.20f, %.20f));\n' % (c_joint.getName(), coord, rot2_axis[0], rot2_axis[1], rot2_axis[2]))
            
            #Rotation 3
            rot3 = spatialtransform.get_rotation3()
            rot3_axis = rot3.get_axis().to_numpy()
            rot3_f = rot3.get_function()
            coord = 2
            if rot3_f.getConcreteClassName() == 'LinearFunction': 
                rot3_f_obj = osim.LinearFunction.safeDownCast(rot3_f)
                rot3_f_slope = rot3_f_obj.getSlope()
                rot3_f_intercept = rot3_f_obj.getIntercept()                
                c_coord = c_joint.get_coordinates(coord)
                c_coord_name = c_coord.getName()
                f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                f.write('\tst_%s[%i].setFunction(new LinearFunction(%.4f, %.4f));\n' % (c_joint.getName(), coord, rot3_f_slope, rot3_f_intercept))
            elif rot3_f.getConcreteClassName() == 'PolynomialFunction':
                f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                rot3_f_obj = osim.PolynomialFunction.safeDownCast(rot3_f)                
                rot3_f_coeffs = rot3_f_obj.getCoefficients().to_numpy()
                c_nCoeffs = rot3_f_coeffs.shape[0]                
                if c_nCoeffs == 2:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot3_f_coeffs[0], rot3_f_coeffs[1]))
                elif c_nCoeffs == 3:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot3_f_coeffs[0], rot3_f_coeffs[1], rot3_f_coeffs[2]))
                elif c_nCoeffs == 4:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot3_f_coeffs[0], rot3_f_coeffs[1], rot3_f_coeffs[2], rot3_f_coeffs[3]))  
                elif c_nCoeffs == 5:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot3_f_coeffs[0], rot3_f_coeffs[1], rot3_f_coeffs[2], rot3_f_coeffs[3], rot3_f_coeffs[4]))                    
                else:
                    raise ValueError("TODO")
                f.write('\tVector st_%s_%i_coeffs_vec(%i); \n' % (c_joint.getName(), coord, c_nCoeffs))
                f.write('\tfor (int i = 0; i < %i; ++i) st_%s_%i_coeffs_vec[i] = st_%s_%i_coeffs[i]; \n' % (c_nCoeffs, c_joint.getName(), coord, c_joint.getName(), coord))
                f.write('\tst_%s[%i].setFunction(new PolynomialFunction(st_%s_%i_coeffs_vec));\n' % (c_joint.getName(), coord, c_joint.getName(), coord))
            elif rot3_f.getConcreteClassName() == 'MultiplierFunction':
                rot3_f_obj = osim.MultiplierFunction.safeDownCast(rot3_f)
                rot3_f_obj_scale = rot3_f_obj.getScale()
                rot3_f_obj_f = rot3_f_obj.getFunction()
                rot3_f_obj_f_name = rot3_f_obj_f.getConcreteClassName()
                if rot3_f_obj_f_name == 'Constant':
                    rot3_f_obj_f_obj = osim.Constant.safeDownCast(rot3_f_obj_f)
                    rot3_f_obj_f_obj_value = rot3_f_obj_f_obj.getValue()
                    f.write('\tst_%s[%i].setFunction(new MultiplierFunction(new Constant(%.20f), %.20f));\n' % (c_joint.getName(), coord, rot3_f_obj_f_obj_value, rot3_f_obj_scale))
                elif rot3_f_obj_f_name == 'PolynomialFunction':
                    f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                    rot3_f_obj_f_obj = osim.PolynomialFunction.safeDownCast(rot3_f_obj_f)
                    rot3_f_obj_f_coeffs = rot3_f_obj_f_obj.getCoefficients().to_numpy()
                    c_nCoeffs = rot3_f_obj_f_coeffs.shape[0]
                    if c_nCoeffs == 2:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot3_f_obj_f_coeffs[0], rot3_f_obj_f_coeffs[1]))
                    elif c_nCoeffs == 3:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot3_f_obj_f_coeffs[0], rot3_f_obj_f_coeffs[1], rot3_f_obj_f_coeffs[2]))
                    elif c_nCoeffs == 4:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot3_f_obj_f_coeffs[0], rot3_f_obj_f_coeffs[1], rot3_f_obj_f_coeffs[2], rot3_f_obj_f_coeffs[3]))  
                    elif c_nCoeffs == 5:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, rot3_f_obj_f_coeffs[0], rot3_f_obj_f_coeffs[1], rot3_f_obj_f_coeffs[2], rot3_f_obj_f_coeffs[3], rot3_f_obj_f_coeffs[4]))                    
                    else:
                        raise ValueError("TODO")
                    f.write('\tVector st_%s_%i_coeffs_vec(%i); \n' % (c_joint.getName(), coord, c_nCoeffs))
                    f.write('\tfor (int i = 0; i < %i; ++i) st_%s_%i_coeffs_vec[i] = st_%s_%i_coeffs[i]; \n' % (c_nCoeffs, c_joint.getName(), coord, c_joint.getName(), coord))
                    f.write('\tst_%s[%i].setFunction(new MultiplierFunction(new PolynomialFunction(st_%s_%i_coeffs_vec), %.20f));\n' % (c_joint.getName(), coord, c_joint.getName(), coord, rot3_f_obj_scale))
                else:
                    raise ValueError("Not supported")
            elif rot3_f.getConcreteClassName() == 'Constant':
                rot3_f_obj = osim.Constant.safeDownCast(rot3_f)
                rot3_f_obj_value = rot3_f_obj.getValue()
                f.write('\tst_%s[%i].setFunction(new Constant(%.20f));\n' % (c_joint.getName(), coord, rot3_f_obj_value))
            else:
                raise ValueError("Not supported")
            f.write('\tst_%s[%i].setAxis(Vec3(%.20f, %.20f, %.20f));\n' % (c_joint.getName(), coord, rot3_axis[0], rot3_axis[1], rot3_axis[2]))
            
            # Translation 1
            tr1 = spatialtransform.get_translation1()
            tr1_axis = tr1.get_axis().to_numpy()
            tr1_f = tr1.get_function()
            coord = 3
            if tr1_f.getConcreteClassName() == 'LinearFunction':    
                tr1_f_obj = osim.LinearFunction.safeDownCast(tr1_f)
                tr1_f_slope = tr1_f_obj.getSlope()
                tr1_f_intercept = tr1_f_obj.getIntercept()                
                c_coord = c_joint.get_coordinates(coord)
                c_coord_name = c_coord.getName()
                f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                f.write('\tst_%s[%i].setFunction(new LinearFunction(%.4f, %.4f));\n' % (c_joint.getName(), coord, tr1_f_slope, tr1_f_intercept))
            elif tr1_f.getConcreteClassName() == 'PolynomialFunction':
                f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                tr1_f_obj = osim.PolynomialFunction.safeDownCast(tr1_f)                
                tr1_f_coeffs = tr1_f_obj.getCoefficients().to_numpy()
                c_nCoeffs = tr1_f_coeffs.shape[0]                
                if c_nCoeffs == 2:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr1_f_coeffs[0], tr1_f_coeffs[1]))
                elif c_nCoeffs == 3:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr1_f_coeffs[0], tr1_f_coeffs[1], tr1_f_coeffs[2]))
                elif c_nCoeffs == 4:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr1_f_coeffs[0], tr1_f_coeffs[1], tr1_f_coeffs[2], tr1_f_coeffs[3]))  
                elif c_nCoeffs == 5:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr1_f_coeffs[0], tr1_f_coeffs[1], tr1_f_coeffs[2], tr1_f_coeffs[3], tr1_f_coeffs[4]))                    
                else:
                    raise ValueError("TODO")
                f.write('\tVector st_%s_%i_coeffs_vec(%i); \n' % (c_joint.getName(), coord, c_nCoeffs))
                f.write('\tfor (int i = 0; i < %i; ++i) st_%s_%i_coeffs_vec[i] = st_%s_%i_coeffs[i]; \n' % (c_nCoeffs, c_joint.getName(), coord, c_joint.getName(), coord))
                f.write('\tst_%s[%i].setFunction(new PolynomialFunction(st_%s_%i_coeffs_vec));\n' % (c_joint.getName(), coord, c_joint.getName(), coord))
            elif tr1_f.getConcreteClassName() == 'MultiplierFunction':
                tr1_f_obj = osim.MultiplierFunction.safeDownCast(tr1_f)
                tr1_f_obj_scale = tr1_f_obj.getScale()
                tr1_f_obj_f = tr1_f_obj.getFunction()
                tr1_f_obj_f_name = tr1_f_obj_f.getConcreteClassName()
                if tr1_f_obj_f_name == 'Constant':
                    tr1_f_obj_f_obj = osim.Constant.safeDownCast(tr1_f_obj_f)
                    tr1_f_obj_f_obj_value = tr1_f_obj_f_obj.getValue()
                    f.write('\tst_%s[%i].setFunction(new MultiplierFunction(new Constant(%.20f), %.20f));\n' % (c_joint.getName(), coord, tr1_f_obj_f_obj_value, tr1_f_obj_scale))
                elif tr1_f_obj_f_name == 'PolynomialFunction':
                    f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                    tr1_f_obj_f_obj = osim.PolynomialFunction.safeDownCast(tr1_f_obj_f)
                    tr1_f_obj_f_coeffs = tr1_f_obj_f_obj.getCoefficients().to_numpy()
                    c_nCoeffs = tr1_f_obj_f_coeffs.shape[0]
                    if c_nCoeffs == 2:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr1_f_obj_f_coeffs[0], tr1_f_obj_f_coeffs[1]))
                    elif c_nCoeffs == 3:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr1_f_obj_f_coeffs[0], tr1_f_obj_f_coeffs[1], tr1_f_obj_f_coeffs[2]))
                    elif c_nCoeffs == 4:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr1_f_obj_f_coeffs[0], tr1_f_obj_f_coeffs[1], tr1_f_obj_f_coeffs[2], tr1_f_obj_f_coeffs[3]))  
                    elif c_nCoeffs == 5:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr1_f_obj_f_coeffs[0], tr1_f_obj_f_coeffs[1], tr1_f_obj_f_coeffs[2], tr1_f_obj_f_coeffs[3], tr1_f_obj_f_coeffs[4]))                    
                    else:
                        raise ValueError("TODO")
                    f.write('\tVector st_%s_%i_coeffs_vec(%i); \n' % (c_joint.getName(), coord, c_nCoeffs))
                    f.write('\tfor (int i = 0; i < %i; ++i) st_%s_%i_coeffs_vec[i] = st_%s_%i_coeffs[i]; \n' % (c_nCoeffs, c_joint.getName(), coord, c_joint.getName(), coord))
                    f.write('\tst_%s[%i].setFunction(new MultiplierFunction(new PolynomialFunction(st_%s_%i_coeffs_vec), %.20f));\n' % (c_joint.getName(), coord, c_joint.getName(), coord, tr1_f_obj_scale))
                else:
                    raise ValueError("Not supported")
            elif tr1_f.getConcreteClassName() == 'Constant':
                tr1_f_obj = osim.Constant.safeDownCast(tr1_f)
                tr1_f_obj_value = tr1_f_obj.getValue()
                f.write('\tst_%s[%i].setFunction(new Constant(%.20f));\n' % (c_joint.getName(), coord, tr1_f_obj_value))
            else:
                raise ValueError("Not supported")
            f.write('\tst_%s[%i].setAxis(Vec3(%.20f, %.20f, %.20f));\n' % (c_joint.getName(), coord, tr1_axis[0], tr1_axis[1], tr1_axis[2]))            
            
            #Translation 2
            tr2 = spatialtransform.get_translation2()
            tr2_axis = tr2.get_axis().to_numpy()
            tr2_f = tr2.get_function()
            coord = 4
            if tr2_f.getConcreteClassName() == 'LinearFunction': 
                tr2_f_obj = osim.LinearFunction.safeDownCast(tr2_f)
                tr2_f_slope = tr2_f_obj.getSlope()
                tr2_f_intercept = tr2_f_obj.getIntercept()                
                c_coord = c_joint.get_coordinates(coord)
                c_coord_name = c_coord.getName()
                f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                f.write('\tst_%s[%i].setFunction(new LinearFunction(%.4f, %.4f));\n' % (c_joint.getName(), coord, tr2_f_slope, tr2_f_intercept))
            elif tr2_f.getConcreteClassName() == 'PolynomialFunction':
                f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                tr2_f_obj = osim.PolynomialFunction.safeDownCast(tr2_f)                
                tr2_f_coeffs = tr2_f_obj.getCoefficients().to_numpy()
                c_nCoeffs = tr2_f_coeffs.shape[0]                
                if c_nCoeffs == 2:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr2_f_coeffs[0], tr2_f_coeffs[1]))
                elif c_nCoeffs == 3:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr2_f_coeffs[0], tr2_f_coeffs[1], tr2_f_coeffs[2]))
                elif c_nCoeffs == 4:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr2_f_coeffs[0], tr2_f_coeffs[1], tr2_f_coeffs[2], tr2_f_coeffs[3]))  
                elif c_nCoeffs == 5:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr2_f_coeffs[0], tr2_f_coeffs[1], tr2_f_coeffs[2], tr2_f_coeffs[3], tr2_f_coeffs[4]))                    
                else:
                    raise ValueError("TODO")
                f.write('\tVector st_%s_%i_coeffs_vec(%i); \n' % (c_joint.getName(), coord, c_nCoeffs))
                f.write('\tfor (int i = 0; i < %i; ++i) st_%s_%i_coeffs_vec[i] = st_%s_%i_coeffs[i]; \n' % (c_nCoeffs, c_joint.getName(), coord, c_joint.getName(), coord))
                f.write('\tst_%s[%i].setFunction(new PolynomialFunction(st_%s_%i_coeffs_vec));\n' % (c_joint.getName(), coord, c_joint.getName(), coord))
            elif tr2_f.getConcreteClassName() == 'MultiplierFunction':
                tr2_f_obj = osim.MultiplierFunction.safeDownCast(tr2_f)
                tr2_f_obj_scale = tr2_f_obj.getScale()
                tr2_f_obj_f = tr2_f_obj.getFunction()
                tr2_f_obj_f_name = tr2_f_obj_f.getConcreteClassName()
                if tr2_f_obj_f_name == 'Constant':
                    tr2_f_obj_f_obj = osim.Constant.safeDownCast(tr2_f_obj_f)
                    tr2_f_obj_f_obj_value = tr2_f_obj_f_obj.getValue()
                    f.write('\tst_%s[%i].setFunction(new MultiplierFunction(new Constant(%.20f), %.20f));\n' % (c_joint.getName(), coord, tr2_f_obj_f_obj_value, tr2_f_obj_scale))
                elif tr2_f_obj_f_name == 'PolynomialFunction':
                    f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                    tr2_f_obj_f_obj = osim.PolynomialFunction.safeDownCast(tr2_f_obj_f)
                    tr2_f_obj_f_coeffs = tr2_f_obj_f_obj.getCoefficients().to_numpy()
                    c_nCoeffs = tr2_f_obj_f_coeffs.shape[0]
                    if c_nCoeffs == 2:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr2_f_obj_f_coeffs[0], tr2_f_obj_f_coeffs[1]))
                    elif c_nCoeffs == 3:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr2_f_obj_f_coeffs[0], tr2_f_obj_f_coeffs[1], tr2_f_obj_f_coeffs[2]))
                    elif c_nCoeffs == 4:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr2_f_obj_f_coeffs[0], tr2_f_obj_f_coeffs[1], tr2_f_obj_f_coeffs[2], tr2_f_obj_f_coeffs[3]))  
                    elif c_nCoeffs == 5:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr2_f_obj_f_coeffs[0], tr2_f_obj_f_coeffs[1], tr2_f_obj_f_coeffs[2], tr2_f_obj_f_coeffs[3], tr2_f_obj_f_coeffs[4]))                    
                    else:
                        raise ValueError("TODO")
                    f.write('\tVector st_%s_%i_coeffs_vec(%i); \n' % (c_joint.getName(), coord, c_nCoeffs))
                    f.write('\tfor (int i = 0; i < %i; ++i) st_%s_%i_coeffs_vec[i] = st_%s_%i_coeffs[i]; \n' % (c_nCoeffs, c_joint.getName(), coord, c_joint.getName(), coord))
                    f.write('\tst_%s[%i].setFunction(new MultiplierFunction(new PolynomialFunction(st_%s_%i_coeffs_vec), %.20f));\n' % (c_joint.getName(), coord, c_joint.getName(), coord, tr2_f_obj_scale))
                else:
                    raise ValueError("Not supported")
            elif tr2_f.getConcreteClassName() == 'Constant':
                tr2_f_obj = osim.Constant.safeDownCast(tr2_f)
                tr2_f_obj_value = tr2_f_obj.getValue()
                f.write('\tst_%s[%i].setFunction(new Constant(%.20f));\n' % (c_joint.getName(), coord, tr2_f_obj_value))
            else:
                raise ValueError("Not supported")
            f.write('\tst_%s[%i].setAxis(Vec3(%.20f, %.20f, %.20f));\n' % (c_joint.getName(), coord, tr2_axis[0], tr2_axis[1], tr2_axis[2]))
            
            # Translation 3
            tr3 = spatialtransform.get_translation3()
            tr3_axis = tr3.get_axis().to_numpy()
            tr3_f = tr3.get_function()
            coord = 5
            if tr3_f.getConcreteClassName() == 'LinearFunction':     
                tr3_f_obj = osim.LinearFunction.safeDownCast(tr3_f)
                tr3_f_slope = tr3_f_obj.getSlope()
                tr3_f_intercept = tr3_f_obj.getIntercept()                
                c_coord = c_joint.get_coordinates(coord)
                c_coord_name = c_coord.getName()
                f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                f.write('\tst_%s[%i].setFunction(new LinearFunction(%.4f, %.4f));\n' % (c_joint.getName(), coord, tr3_f_slope, tr3_f_intercept))
            elif tr3_f.getConcreteClassName() == 'PolynomialFunction':
                f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                tr3_f_obj = osim.PolynomialFunction.safeDownCast(tr3_f)                
                tr3_f_coeffs = tr3_f_obj.getCoefficients().to_numpy()
                c_nCoeffs = tr3_f_coeffs.shape[0]                
                if c_nCoeffs == 2:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr3_f_coeffs[0], tr3_f_coeffs[1]))
                elif c_nCoeffs == 3:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr3_f_coeffs[0], tr3_f_coeffs[1], tr3_f_coeffs[2]))
                elif c_nCoeffs == 4:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr3_f_coeffs[0], tr3_f_coeffs[1], tr3_f_coeffs[2], tr3_f_coeffs[3]))  
                elif c_nCoeffs == 5:
                    f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr3_f_coeffs[0], tr3_f_coeffs[1], tr3_f_coeffs[2], tr3_f_coeffs[3], tr3_f_coeffs[4]))                    
                else:
                    raise ValueError("TODO")
                f.write('\tVector st_%s_%i_coeffs_vec(%i); \n' % (c_joint.getName(), coord, c_nCoeffs))
                f.write('\tfor (int i = 0; i < %i; ++i) st_%s_%i_coeffs_vec[i] = st_%s_%i_coeffs[i]; \n' % (c_nCoeffs, c_joint.getName(), coord, c_joint.getName(), coord))
                f.write('\tst_%s[%i].setFunction(new PolynomialFunction(st_%s_%i_coeffs_vec));\n' % (c_joint.getName(), coord, c_joint.getName(), coord))
            elif tr3_f.getConcreteClassName() == 'MultiplierFunction':
                tr3_f_obj = osim.MultiplierFunction.safeDownCast(tr3_f)
                tr3_f_obj_scale = tr3_f_obj.getScale()
                tr3_f_obj_f = tr3_f_obj.getFunction()
                tr3_f_obj_f_name = tr3_f_obj_f.getConcreteClassName()
                if tr3_f_obj_f_name == 'Constant':
                    tr3_f_obj_f_obj = osim.Constant.safeDownCast(tr3_f_obj_f)
                    tr3_f_obj_f_obj_value = tr3_f_obj_f_obj.getValue()
                    f.write('\tst_%s[%i].setFunction(new MultiplierFunction(new Constant(%.20f), %.20f));\n' % (c_joint.getName(), coord, tr3_f_obj_f_obj_value, tr3_f_obj_scale))
                elif tr3_f_obj_f_name == 'PolynomialFunction':
                    f.write('\tst_%s[%i].setCoordinateNames(OpenSim::Array<std::string>(\"%s\", 1, 1));\n' % (c_joint.getName(), coord, c_coord_name))
                    tr3_f_obj_f_obj = osim.PolynomialFunction.safeDownCast(tr3_f_obj_f)
                    tr3_f_obj_f_coeffs = tr3_f_obj_f_obj.getCoefficients().to_numpy()
                    c_nCoeffs = tr3_f_obj_f_coeffs.shape[0]
                    if c_nCoeffs == 2:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr3_f_obj_f_coeffs[0], tr3_f_obj_f_coeffs[1]))
                    elif c_nCoeffs == 3:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr3_f_obj_f_coeffs[0], tr3_f_obj_f_coeffs[1], tr3_f_obj_f_coeffs[2]))
                    elif c_nCoeffs == 4:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr3_f_obj_f_coeffs[0], tr3_f_obj_f_coeffs[1], tr3_f_obj_f_coeffs[2], tr3_f_obj_f_coeffs[3]))  
                    elif c_nCoeffs == 5:
                        f.write('\tosim_double_adouble st_%s_%i_coeffs[%i] = {%.20f, %.20f, %.20f, %.20f, %.20f}; \n' % (c_joint.getName(), coord, c_nCoeffs, tr3_f_obj_f_coeffs[0], tr3_f_obj_f_coeffs[1], tr3_f_obj_f_coeffs[2], tr3_f_obj_f_coeffs[3], tr3_f_obj_f_coeffs[4]))                    
                    else:
                        raise ValueError("TODO")
                    f.write('\tVector st_%s_%i_coeffs_vec(%i); \n' % (c_joint.getName(), coord, c_nCoeffs))
                    f.write('\tfor (int i = 0; i < %i; ++i) st_%s_%i_coeffs_vec[i] = st_%s_%i_coeffs[i]; \n' % (c_nCoeffs, c_joint.getName(), coord, c_joint.getName(), coord))
                    f.write('\tst_%s[%i].setFunction(new MultiplierFunction(new PolynomialFunction(st_%s_%i_coeffs_vec), %.20f));\n' % (c_joint.getName(), coord, c_joint.getName(), coord, tr3_f_obj_scale))
                else:
                    raise ValueError("Not supported") 
            elif tr3_f.getConcreteClassName() == 'Constant':
                tr3_f_obj = osim.Constant.safeDownCast(tr3_f)
                tr3_f_obj_value = tr3_f_obj.getValue()
                f.write('\tst_%s[%i].setFunction(new Constant(%.20f));\n' % (c_joint.getName(), coord, tr3_f_obj_value))
            else:
                raise ValueError("Not supported")
            f.write('\tst_%s[%i].setAxis(Vec3(%.20f, %.20f, %.20f));\n' % (c_joint.getName(), coord, tr3_axis[0], tr3_axis[1], tr3_axis[2]))          
            
            #Joint.
            f.write('\tOpenSim::%s* %s;\n' % (c_joint_type, c_joint.getName()))
            if parent_frame_name == "ground":
                f.write('\t%s = new OpenSim::%s(\"%s\", model->getGround(), Vec3(%.20f, %.20f, %.20f), Vec3(%.20f, %.20f, %.20f), *%s, Vec3(%.20f, %.20f, %.20f), Vec3(%.20f, %.20f, %.20f), st_%s);\n' % (c_joint.getName(), c_joint_type, c_joint.getName(), parent_frame_trans[0], parent_frame_trans[1], parent_frame_trans[2], parent_frame_or[0], parent_frame_or[1], parent_frame_or[2], child_frame_name, child_frame_trans[0], child_frame_trans[1], child_frame_trans[2], child_frame_or[0], child_frame_or[1], child_frame_or[2], c_joint.getName()))     
            else:
                f.write('\t%s = new OpenSim::%s(\"%s\", *%s, Vec3(%.20f, %.20f, %.20f), Vec3(%.20f, %.20f, %.20f), *%s, Vec3(%.20f, %.20f, %.20f), Vec3(%.20f, %.20f, %.20f), st_%s);\n' % (c_joint.getName(), c_joint_type, c_joint.getName(), parent_frame_name, parent_frame_trans[0], parent_frame_trans[1], parent_frame_trans[2], parent_frame_or[0], parent_frame_or[1], parent_frame_or[2], child_frame_name, child_frame_trans[0], child_frame_trans[1], child_frame_trans[2], child_frame_or[0], child_frame_or[1], child_frame_or[2], c_joint.getName()))
            
        elif c_joint_type == 'PinJoint' or c_joint_type == 'WeldJoint' :
            f.write('\tOpenSim::%s* %s;\n' % (c_joint_type, c_joint.getName()))
            if parent_frame_name == "ground":
                f.write('\t%s = new OpenSim::%s(\"%s\", model->getGround(), Vec3(%.20f, %.20f, %.20f), Vec3(%.20f, %.20f, %.20f), *%s, Vec3(%.20f, %.20f, %.20f), Vec3(%.20f, %.20f, %.20f));\n' % (c_joint.getName(), c_joint_type, c_joint.getName(), parent_frame_trans[0], parent_frame_trans[1], parent_frame_trans[2], parent_frame_or[0], parent_frame_or[1], parent_frame_or[2], child_frame_name, child_frame_trans[0], child_frame_trans[1], child_frame_trans[2], child_frame_or[0], child_frame_or[1], child_frame_or[2]))     
            else:
                f.write('\t%s = new OpenSim::%s(\"%s\", *%s, Vec3(%.20f, %.20f, %.20f), Vec3(%.20f, %.20f, %.20f), *%s, Vec3(%.20f, %.20f, %.20f), Vec3(%.20f, %.20f, %.20f));\n' % (c_joint.getName(), c_joint_type, c_joint.getName(), parent_frame_name, parent_frame_trans[0], parent_frame_trans[1], parent_frame_trans[2], parent_frame_or[0], parent_frame_or[1], parent_frame_or[2], child_frame_name, child_frame_trans[0], child_frame_trans[1], child_frame_trans[2], child_frame_or[0], child_frame_or[1], child_frame_or[2])) 
        else:
            raise ValueError("TODO: joint type not yet supported")
        f.write('\tmodel->addJoint(%s);\n' % (c_joint.getName()))
        f.write('\n')  
    if treadmill:
        f.write('\tOpenSim::SliderJoint* ground_treadmill;\n')
        f.write('\tground_treadmill = new SliderJoint("ground_treadmill", model->getGround(), Vec3(0), Vec3(0), *treadmill, Vec3(0), Vec3(0));\n')
        f.write('\tmodel->addJoint(ground_treadmill);\n')
        f.write('\n')
        
    #Contacts
    f.write('\t// Definition of contacts.\n')
    rightFootContact = False
    leftFootContact = False
    rightFootContactBodies = []
    leftFootContactBodies = []
    nRightContacts = 0
    nLeftContacts = 0
    for i in range(forceSet.getSize()):        
        c_force_elt = forceSet.get(i)
        if c_force_elt.getConcreteClassName() == "SmoothSphereHalfSpaceForce":            
            c_force_elt_obj =  osim.SmoothSphereHalfSpaceForce.safeDownCast(c_force_elt) 	
            
            socket0Name = c_force_elt.getSocketNames()[0]
            socket0 = c_force_elt.getSocket(socket0Name)
            socket0_obj = socket0.getConnecteeAsObject()
            socket0_objName = socket0_obj.getName()            
            geo0 = geometrySet.get(socket0_objName)
            geo0_loc = geo0.get_location().to_numpy()
            geo0_or = geo0.get_orientation().to_numpy()
            geo0_frameName = geo0.getFrame().getName()
            
            socket1Name = c_force_elt.getSocketNames()[1]
            socket1 = c_force_elt.getSocket(socket1Name)
            socket1_obj = socket1.getConnecteeAsObject()
            socket1_objName = socket1_obj.getName()            
            geo1 = geometrySet.get(socket1_objName)
            geo1_loc = geo1.get_location().to_numpy()
            geo1_frameName = geo1.getFrame().getName()
            obj = osim.ContactSphere.safeDownCast(geo1) 	
            geo1_radius = obj.getRadius()            
            
            f.write('\tOpenSim::%s* %s;\n' % (c_force_elt.getConcreteClassName(), c_force_elt.getName()))
            if geo0_frameName == "ground":
                if treadmill:
                    ground_contact = "*treadmill"
                else:
                    ground_contact = "model->getGround()"
                
                f.write('\t%s = new %s(\"%s\", *%s, %s);\n' % (c_force_elt.getName(), c_force_elt.getConcreteClassName(), c_force_elt.getName(), geo1_frameName, ground_contact))
            else:
                f.write('\t%s = new %s(\"%s\", *%s, *%s);\n' % (c_force_elt.getName(), c_force_elt.getConcreteClassName(), c_force_elt.getName(), geo1_frameName, geo0_frameName))
                
            f.write('\tVec3 %s_location(%.20f, %.20f, %.20f);\n' % (c_force_elt.getName(), geo1_loc[0], geo1_loc[1], geo1_loc[2]))
            f.write('\t%s->set_contact_sphere_location(%s_location);\n' % (c_force_elt.getName(), c_force_elt.getName()))
            f.write('\tdouble %s_radius = (%.20f);\n' % (c_force_elt.getName(), geo1_radius))
            f.write('\t%s->set_contact_sphere_radius(%s_radius );\n' % (c_force_elt.getName(), c_force_elt.getName()))
            f.write('\t%s->set_contact_half_space_location(Vec3(%.20f, %.20f, %.20f));\n' % (c_force_elt.getName(), geo0_loc[0], geo0_loc[1], geo0_loc[2]))
            f.write('\t%s->set_contact_half_space_orientation(Vec3(%.20f, %.20f, %.20f));\n' % (c_force_elt.getName(), geo0_or[0], geo0_or[1], geo0_or[2]))
            
            f.write('\t%s->set_stiffness(%.20f);\n' % (c_force_elt.getName(), c_force_elt_obj.get_stiffness()))
            f.write('\t%s->set_dissipation(%.20f);\n' % (c_force_elt.getName(), c_force_elt_obj.get_dissipation()))
            f.write('\t%s->set_static_friction(%.20f);\n' % (c_force_elt.getName(), c_force_elt_obj.get_static_friction()))
            f.write('\t%s->set_dynamic_friction(%.20f);\n' % (c_force_elt.getName(), c_force_elt_obj.get_dynamic_friction()))
            f.write('\t%s->set_viscous_friction(%.20f);\n' % (c_force_elt.getName(), c_force_elt_obj.get_viscous_friction()))
            f.write('\t%s->set_transition_velocity(%.20f);\n' % (c_force_elt.getName(), c_force_elt_obj.get_transition_velocity()))
            
            f.write('\t%s->connectSocket_sphere_frame(*%s);\n' % (c_force_elt.getName(), geo1_frameName))
            if geo0_frameName == "ground":
                f.write('\t%s->connectSocket_half_space_frame(%s);\n' % (c_force_elt.getName(), ground_contact))                
            else:
                f.write('\t%s->connectSocket_half_space_frame(*%s);\n' % (c_force_elt.getName(), geo0_frameName))
            f.write('\tmodel->addComponent(%s);\n' % (c_force_elt.getName()))
            f.write('\n')

            #Check if there are right and left foot contacts
            if c_force_elt.getName()[-2:] == '_r':
                nRightContacts += 1
                rightFootContactBodies.append(geo1_frameName)
                if not rightFootContact:
                    rightFootContact = True
            if c_force_elt.getName()[-2:] == '_l':
                nLeftContacts += 1
                leftFootContactBodies.append(geo1_frameName)
                if not leftFootContact:
                    leftFootContact = True
    nContacts = nRightContacts + nLeftContacts
       
    #Compute residuals (joint torques).
    f.write('\t// Initialize system.\n')
    f.write('\tSimTK::State* state;\n')
    f.write('\tstate = new State(model->initSystem());\n\n')

    f.write('\t// Read inputs.\n')
    f.write('\tstd::vector<T> x(arg[0], arg[0] + NX);\n')
    f.write('\tstd::vector<T> u(arg[1], arg[1] + NU);\n')
    if treadmill:
        f.write('\tstd::vector<T> p(arg[2], arg[2] + 1);\n')
    f.write('\n')
    
    f.write('\t// States and controls.\n')
    if treadmill:
        f.write('\tT ua[NU_treadmill];\n')
        f.write('\tVector QsUs(NX_treadmill);\n')
    else:
        f.write('\tT ua[NU];\n')
        f.write('\tVector QsUs(NX);\n')        
    f.write('\t/// States\n')
    f.write('\tfor (int i = 0; i < NX; ++i) QsUs[i] = x[i];\n')
    if treadmill:
        f.write('\tQsUs[NX] = 0;\n')
        f.write('\tQsUs[NX+1] = p[0];\n')
    f.write('\t/// Controls\n')
    if treadmill:
        f.write('\tT ut[NU_treadmill];\n')
        f.write('\tfor (int i = 0; i < NU; ++i) ut[i] = u[i];\n')
        f.write('\tut[NU] = 0;\n')        
    f.write('\t/// OpenSim and Simbody have different state orders.\n')
    f.write('\tauto indicesOSInSimbody = getIndicesOpenSimInSimbody(*model);\n')
    if treadmill:
        f.write('\tfor (int i = 0; i < NU_treadmill; ++i) ua[i] = ut[indicesOSInSimbody[i]];\n\n')
    else:
        f.write('\tfor (int i = 0; i < NU; ++i) ua[i] = u[indicesOSInSimbody[i]];\n\n')

    f.write('\t// Set state variables and realize.\n')
    f.write('\tmodel->setStateVariableValues(*state, QsUs);\n')
    f.write('\tmodel->realizeVelocity(*state);\n\n')
    
    f.write('\t// Compute residual forces.\n')
    f.write('\t/// Set appliedMobilityForces (# mobilities).\n')
    if treadmill:
        f.write('\tVector appliedMobilityForces(nCoordinates_treadmill);\n')
    else:
        f.write('\tVector appliedMobilityForces(nCoordinates);\n')
    f.write('\tappliedMobilityForces.setToZero();\n')
    f.write('\t/// Set appliedBodyForces (# bodies + ground).\n')
    f.write('\tVector_<SpatialVec> appliedBodyForces;\n')
    f.write('\tint nbodies = model->getBodySet().getSize() + 1;\n')
    f.write('\tappliedBodyForces.resize(nbodies);\n')
    f.write('\tappliedBodyForces.setToZero();\n')
    f.write('\t/// Set gravity.\n')
    f.write('\tVec3 gravity(0);\n')
    f.write('\tgravity[1] = %.20f;\n' % model.get_gravity()[1])
    f.write('\t/// Add weights to appliedBodyForces.\n')
    f.write('\tfor (int i = 0; i < model->getBodySet().getSize(); ++i) {\n')
    f.write('\t\tmodel->getMatterSubsystem().addInStationForce(*state,\n')
    f.write('\t\tmodel->getBodySet().get(i).getMobilizedBodyIndex(),\n')
    f.write('\t\tmodel->getBodySet().get(i).getMassCenter(),\n')
    f.write('\t\tmodel->getBodySet().get(i).getMass()*gravity, appliedBodyForces);\n')
    f.write('\t}\n')    
    f.write('\t/// Add contact forces to appliedBodyForces.\n')
    
    count = 0
    for i in range(forceSet.getSize()):        
        c_force_elt = forceSet.get(i)     
        
        if c_force_elt.getConcreteClassName() == "SmoothSphereHalfSpaceForce":
            c_force_elt_name = c_force_elt.getName()    
            
            f.write('\tArray<osim_double_adouble> Force_%s = %s->getRecordValues(*state);\n' % (str(count), c_force_elt_name))
            f.write('\tSpatialVec GRF_%s;\n' % (str(count)))           
            
            f.write('\tGRF_%s[0] = Vec3(Force_%s[3], Force_%s[4], Force_%s[5]);\n' % (str(count), str(count), str(count), str(count)))
            f.write('\tGRF_%s[1] = Vec3(Force_%s[0], Force_%s[1], Force_%s[2]);\n' % (str(count), str(count), str(count), str(count)))
            
            socket1Name = c_force_elt.getSocketNames()[1]
            socket1 = c_force_elt.getSocket(socket1Name)
            socket1_obj = socket1.getConnecteeAsObject()
            socket1_objName = socket1_obj.getName()            
            geo1 = geometrySet.get(socket1_objName)
            geo1_frameName = geo1.getFrame().getName()
            
            f.write('\tint c_idx_%s = model->getBodySet().get("%s").getMobilizedBodyIndex();\n' % (str(count), geo1_frameName))            
            f.write('\tappliedBodyForces[c_idx_%s] += GRF_%s;\n' % (str(count), str(count)))
            count += 1
            f.write('\n')
            
    f.write('\t/// knownUdot.\n')
    if treadmill:
        f.write('\tVector knownUdot(nCoordinates_treadmill);\n')
    else:
        f.write('\tVector knownUdot(nCoordinates);\n')
    f.write('\tknownUdot.setToZero();\n')
    if treadmill:
        f.write('\tfor (int i = 0; i < nCoordinates_treadmill; ++i) knownUdot[i] = ua[i];\n')
    else:
        f.write('\tfor (int i = 0; i < nCoordinates; ++i) knownUdot[i] = ua[i];\n')
    f.write('\t/// Calculate residual forces.\n')
    if treadmill:
        f.write('\tVector residualMobilityForces(nCoordinates_treadmill);\n')
    else:
        f.write('\tVector residualMobilityForces(nCoordinates);\n')
    f.write('\tresidualMobilityForces.setToZero();\n')
    f.write('\tmodel->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,\n')
    f.write('\t\t\tappliedMobilityForces, appliedBodyForces, knownUdot, residualMobilityForces);\n\n')
        
    #Get body origins.
    f.write('\t/// Body origins.\n')
    for i in range(bodySet.getSize()):        
        c_body = bodySet.get(i)
        c_body_name = c_body.getName()            
        if (c_body_name == 'patella_l' or c_body_name == 'patella_r'):
            continue            
        f.write('\tVec3 %s_or = %s->getPositionInGround(*state);\n' % (c_body_name, c_body_name))
    f.write('\n')
        
    #Get GRFs.
    f.write('\t/// Ground reaction forces.\n')
    if rightFootContact:
        f.write('\tVec3 GRF_r(0);\n')
    if leftFootContact:
        f.write('\tVec3 GRF_l(0);\n')
    count = 0
    for i in range(forceSet.getSize()):        
        c_force_elt = forceSet.get(i)  
        if c_force_elt.getConcreteClassName() == "SmoothSphereHalfSpaceForce":
            c_force_elt_name = c_force_elt.getName() 
            if c_force_elt_name[-2:] == "_r":
                f.write('\tGRF_r += GRF_%s[1];\n'  % (str(count)))
            elif c_force_elt_name[-2:] == "_l":
                f.write('\tGRF_l += GRF_%s[1];\n'  % (str(count)))
            else:
                raise ValueError("Cannot identify contact side")
            count += 1
    f.write('\n')
        
    #Get GRMs.
    f.write('\t/// Ground reaction moments.\n')
    if rightFootContact:
        f.write('\tVec3 GRM_r(0);\n')
    if leftFootContact:
        f.write('\tVec3 GRM_l(0);\n')
    f.write('\tVec3 normal(0, 1, 0);\n\n')
    count = 0
    geo1_frameNames = []
    for i in range(forceSet.getSize()):        
        c_force_elt = forceSet.get(i)  
        if c_force_elt.getConcreteClassName() == "SmoothSphereHalfSpaceForce":
            c_force_elt_name = c_force_elt.getName() 
            socket1Name = c_force_elt.getSocketNames()[1]
            socket1 = c_force_elt.getSocket(socket1Name)
            socket1_obj = socket1.getConnecteeAsObject()
            socket1_objName = socket1_obj.getName()            
            geo1 = geometrySet.get(socket1_objName)
            geo1_frameName = geo1.getFrame().getName() 
            
            if not geo1_frameName in geo1_frameNames:
                f.write('\tSimTK::Transform TR_GB_%s = %s->getMobilizedBody().getBodyTransform(*state);\n' % (geo1_frameName, geo1_frameName))    
                geo1_frameNames.append(geo1_frameName)
                
            f.write('\tVec3 %s_location_G = %s->findStationLocationInGround(*state, %s_location);\n' % (c_force_elt_name, geo1_frameName, c_force_elt_name))                
            f.write('\tVec3 %s_locationCP_G = %s_location_G - %s_radius * normal;\n' % (c_force_elt_name, c_force_elt_name, c_force_elt_name))
            f.write('\tVec3 locationCP_G_adj_%i = %s_locationCP_G - 0.5*%s_locationCP_G[1] * normal;\n' % (count, c_force_elt_name, c_force_elt_name))
            f.write('\tVec3 %s_locationCP_B = model->getGround().findStationLocationInAnotherFrame(*state, locationCP_G_adj_%i, *%s);\n' % (c_force_elt_name, count, geo1_frameName))
            f.write('\tVec3 GRM_%i = (TR_GB_%s*%s_locationCP_B) %% GRF_%s[1];\n' % (count, geo1_frameName, c_force_elt_name, str(count)))
            
            if c_force_elt_name[-2:] == "_r":
                f.write('\tGRM_r += GRM_%i;\n'  % (count))   
            elif c_force_elt_name[-2:] == "_l": 
                f.write('\tGRM_l += GRM_%i;\n'  % (count))   
            else:
                raise ValueError("Cannot identify contact side")
            f.write('\n')                   
            count += 1
    
    #Save dict pointing to which elements are returned by F and in which
    #order, such as to facilitate using F when formulating problem.
    F_map = {}
    
    f.write('\t/// Outputs.\n')        
    
    #Export residuals (joint torques).
    f.write('\t/// Residual forces (OpenSim and Simbody have different state orders).\n')
    f.write('\tauto indicesSimbodyInOS = getIndicesSimbodyInOpenSim(*model);\n')
    f.write('\tfor (int i = 0; i < NU; ++i) res[0][i] =\n')
    f.write('\t\t\tvalue<T>(residualMobilityForces[indicesSimbodyInOS[i]]);\n')
    F_map['residuals'] = {}
    count = 0
    for coordinate in coordinates:
        if 'beta' in coordinate:
            continue
        F_map['residuals'][coordinate] = count 
        count += 1
    count_acc = nCoordinates
    
    #Export GRFs.
    f.write('\t/// Ground reaction forces.\n')        
    F_map['GRFs'] = {} 
    F_map['GRFs']['nContactSpheres'] = nContacts
    F_map['GRFs']['nRightContactSpheres'] = nRightContacts
    F_map['GRFs']['nLeftContactSpheres'] = nLeftContacts
    if rightFootContact:
        f.write('\tfor (int i = 0; i < 3; ++i) res[0][i + %i] = value<T>(GRF_r[i]);\n' % (count_acc))
        F_map['GRFs']['right'] = range(count_acc, count_acc+3)
        count_acc += 3
    if leftFootContact:
        f.write('\tfor (int i = 0; i < 3; ++i) res[0][i + %i] = value<T>(GRF_l[i]);\n' % (count_acc))
        F_map['GRFs']['left'] = range(count_acc, count_acc+3)
        count_acc += 3       
    
    #Export GRMs.
    f.write('\t/// Ground reaction moments.\n')
    F_map['GRMs'] = {}
    if rightFootContact:
        f.write('\tfor (int i = 0; i < 3; ++i) res[0][i + %i] = value<T>(GRM_r[i]);\n' % (count_acc))
        F_map['GRMs']['right'] = range(count_acc, count_acc+3)
        count_acc += 3
    if leftFootContact:
        f.write('\tfor (int i = 0; i < 3; ++i) res[0][i + %i] = value<T>(GRM_l[i]);\n' % (count_acc))
        F_map['GRMs']['left'] = range(count_acc, count_acc+3)
        count_acc += 3
    
    #Export individual GRFs.
    f.write('\t/// Ground reaction forces per sphere.\n')
    count = 0
    F_map['GRFs']['rightContactSpheres'] = []
    F_map['GRFs']['leftContactSpheres'] = []
    F_map['GRFs']['rightContactSphereBodies'] = rightFootContactBodies
    F_map['GRFs']['leftContactSphereBodies'] = leftFootContactBodies        
    for i in range(forceSet.getSize()):
        c_force_elt = forceSet.get(i) 
        if c_force_elt.getConcreteClassName() == "SmoothSphereHalfSpaceForce":
            f.write('\tfor (int i = 0; i < 3; ++i) res[0][i + %i] = value<T>(GRF_%i[1][i]);\n' % (count_acc, count))
            F_map['GRFs'][c_force_elt.getName()] = range(count_acc, count_acc+3)
            if c_force_elt.getName()[-2:] == "_r":
                F_map['GRFs']['rightContactSpheres'].append(c_force_elt.getName())
            elif c_force_elt.getName()[-2:] == "_l":
                F_map['GRFs']['leftContactSpheres'].append(c_force_elt.getName())
            count += 1
            count_acc += 3
    f.write('\n')
    
    #Export individual contact locations.
    f.write('\t/// Contact point locations per sphere.\n')
    F_map['COPs'] = {}
    count = 0
    for i in range(forceSet.getSize()):
        c_force_elt = forceSet.get(i) 
        if c_force_elt.getConcreteClassName() == "SmoothSphereHalfSpaceForce":
            f.write('\tfor (int i = 0; i < 3; ++i) res[0][i + %i] = value<T>(locationCP_G_adj_%i[i]);\n' % (count_acc, count))
            F_map['COPs'][c_force_elt.getName()] = range(count_acc, count_acc+3)
            count += 1
            count_acc += 3
    f.write('\n')
    
    #Export body origins.
    f.write('\t/// Body origins.\n')
    F_map['body_origins'] = {}
    count = 0
    for i in range(bodySet.getSize()):        
        c_body = bodySet.get(i)
        c_body_name = c_body.getName()
        if (c_body_name == 'patella_l' or c_body_name == 'patella_r'):
            continue
        f.write('\tfor (int i = 0; i < 3; ++i) res[0][i + %i] = value<T>(%s_or[i]);\n' % (count_acc+count*3, c_body_name))
        F_map['body_origins'][c_body_name] = range(count_acc+count*3, count_acc+count*3+3)
        count += 1
    count_acc += 3*count
        
    f.write('\n')
    f.write('\treturn 0;\n')
    f.write('}\n\n')
    
    #Residuals (joint torques), 3D GRFs (combined), 3D GRMs (combined),
    #3D GRFs (per sphere), 3D COP (per sphere), and 3D body origins.
    nOutputs = nCoordinates + 3*(2*nContacts+nBodies)
    if rightFootContact:
        nOutputs += 2*3
    if leftFootContact:
        nOutputs += 2*3
    f.write('constexpr int NR = %i; \n\n' % (nOutputs))
    
    f.write('int main() {\n')
    f.write('\tRecorder x[NX];\n')
    f.write('\tRecorder u[NU];\n')
    if treadmill:
        f.write('\tRecorder p[1];\n')            
    f.write('\tRecorder tau[NR];\n')
    f.write('\tfor (int i = 0; i < NX; ++i) x[i] <<= 0;\n')
    f.write('\tfor (int i = 0; i < NU; ++i) u[i] <<= 0;\n')
    if treadmill:
        f.write('\tp[0] <<= 0;\n')
        f.write('\tconst Recorder* Recorder_arg[n_in] = { x,u,p };\n')
    else:
        f.write('\tconst Recorder* Recorder_arg[n_in] = { x,u };\n')
    f.write('\tRecorder* Recorder_res[n_out] = { tau };\n')
    f.write('\tF_generic<Recorder>(Recorder_arg, Recorder_res);\n')
    f.write('\tdouble res[NR];\n')
    f.write('\tfor (int i = 0; i < NR; ++i) Recorder_res[0][i] >>= res[i];\n')
    f.write('\tRecorder::stop_recording();\n')
    f.write('\treturn 0;\n')
    f.write('}\n')
    
    #Save dict.
    np.save(pathOutputMap, F_map)
        
#Build external Function

"""
Note that this build external functions requires the functions included in the
UtilsDynamicSimulations folder

It also takes the buildExternalFunction function from the utilsOpenSimAD script
in the opencap-processing example.

buildExternalFunction(
    externalFunctionName, pathDCAD, pathOutputExternalFunctionFolder,
    3*nCoordinates, treadmill=treadmill)

"""

#Set utils directory to appropriate folder
baseDir = os.getcwd()+'\\supplementary'

#Set to build external function
build_externalFunction = True

#Check and build external function
if build_externalFunction:
    
    #Set path to dynamic simulation utilities
    pathDCAD = os.path.join(baseDir, 'UtilsDynamicSimulations', 'OpenSimAD') 
    
    #Build external function
    
    #Set inputs to below function commands
    filename = externalFunctionName
    CPP_DIR = pathOutputExternalFunctionFolder
    nInputs = 3*nCoordinates
    
    #Part 1: build expression graph (i.e., generate foo.py).
    pathMain = os.getcwd()
    pathBuildExpressionGraph = os.path.join(pathDCAD, 'buildExpressionGraph')
    pathBuild = os.path.join(pathDCAD, 'build-ExpressionGraph' + filename)
    os.makedirs(pathBuild, exist_ok=True)
    OpenSimAD_DIR = os.path.join(pathDCAD, 'opensimAD-install')
    os.makedirs(OpenSimAD_DIR, exist_ok=True)
    os_system = platform.system()
    
    #Note that Linux and Darwin system options removed
    #The below code is just the Windows option
    
    #Set install, bin and sdk directories
    OpenSimADOS_DIR = os.path.join(OpenSimAD_DIR, 'windows')        
    BIN_DIR = os.path.join(OpenSimADOS_DIR, 'bin')
    SDK_DIR = os.path.join(OpenSimADOS_DIR, 'sdk')
    
    #Download libraries if not existing locally.
    if not os.path.exists(BIN_DIR):
        url = 'https://sourceforge.net/projects/opensimad/files/windows.zip'
        zipfilename = 'windows.zip'
        try:
            download_file(url, zipfilename)
        except:
            try:
                download_file_2(url, zipfilename)
            except:
                error_msg = """ \n\n\n
                Problem when downloading third-party libraries. You can download them manually:
                    1. Download the zip file hosted here: {},
                    2. Extract the files, and
                    3. Copy then under: <local_path>/opencap-processing/UtilsDynamicSimulations/OpenSimAD/opensimAD-install.
                You should have:
                    1. <local_path>/opencap-processing/UtilsDynamicSimulations/OpenSimAD/opensimAD-install/windows/bin and
                    2. <local_path>/opencap-processing/UtilsDynamicSimulations/OpenSimAD/opensimAD-install/windows/sdk \n\n\n""".format(url)
                raise ValueError(error_msg)                    
        with zipfile.ZipFile('windows.zip', 'r') as zip_ref:
            zip_ref.extractall(OpenSimAD_DIR)
        os.remove('windows.zip')
        
    #Create command strings
    cmd1 = 'cmake "' + pathBuildExpressionGraph + '"  -A x64 -DTARGET_NAME:STRING="' + filename + '" -DSDK_DIR:PATH="' + SDK_DIR + '" -DCPP_DIR:PATH="' + CPP_DIR + '"'
    cmd2 = "cmake --build . --config RelWithDebInfo"
        
    #Change to build path    
    os.chdir(pathBuild)  
    
    #Run system commands
    os.system(cmd1)
    os.system(cmd2)
    
    #Create Windows executable
    os.chdir(BIN_DIR)
    path_EXE = os.path.join(pathBuild, 'RelWithDebInfo', filename + '.exe')
    cmd2w = '"{}"'.format(path_EXE)
    os.system(cmd2w)
    
    #Part 2: build external function (i.e., build .dll/.so/.dylib).
    
    #File settings
    fooName = "foo.py"
    pathBuildExternalFunction = os.path.join(pathDCAD, 'buildExternalFunction')
    path_external_filename_foo = os.path.join(BIN_DIR, fooName)
    path_external_functions_filename_build = os.path.join(pathDCAD, 'build-ExternalFunction' + filename)
    path_external_functions_filename_install = os.path.join(pathDCAD, 'install-ExternalFunction' + filename)
    os.makedirs(path_external_functions_filename_build, exist_ok=True) 
    os.makedirs(path_external_functions_filename_install, exist_ok=True)
    shutil.copy2(path_external_filename_foo, pathBuildExternalFunction)
    sys.path.append(pathBuildExternalFunction)
    os.chdir(pathBuildExternalFunction)
    
    #Set the dimensions for generating functions    
    if treadmill:
        dim = nInputs+1
    else:
        dim = nInputs
        
    #Generate C code from expression graph
    """
    Note that this small section comes from the generateF function in utilsOpenSimAD
    """
    import foo
    importlib.reload(foo)
    cg = ca.CodeGenerator('foo_jac')
    arg = ca.SX.sym('arg', dim)
    y,_,_ = foo.foo(arg)
    F = ca.Function('F',[arg],[y])
    cg.add(F)
    cg.add(F.jacobian())
    cg.generate()
    
    #Create Windows system commands
    cmd3 = 'cmake "' + pathBuildExternalFunction + '" -A x64 -DTARGET_NAME:STRING="' + filename + '" -DINSTALL_DIR:PATH="' + path_external_functions_filename_install + '"'
    cmd4 = "cmake --build . --config RelWithDebInfo --target install"
    
    #Change to external functions build directory
    os.chdir(path_external_functions_filename_build)
    
    #Run system commands
    os.system(cmd3)
    os.system(cmd4)
    
    #Return to main path
    os.chdir(pathMain)
    
    #Copy dll to main model directory
    shutil.copy2(os.path.join(path_external_functions_filename_install, 'bin', filename + '.dll'), CPP_DIR)
    
    #Clean-up the mess that was made to create the external function    
    os.remove(os.path.join(pathBuildExternalFunction, "foo_jac.c"))
    os.remove(os.path.join(pathBuildExternalFunction, fooName))
    os.remove(path_external_filename_foo)
    shutil.rmtree(pathBuild)
    shutil.rmtree(path_external_functions_filename_install)
    shutil.rmtree(path_external_functions_filename_build)
    
#Set whether to verify with inverse dynamics
verifyID = True

"""

Note that verification uses the OpenSimPipelines directory that comes with the
opencap-processing example and is in the supplementary folder here.

"""
    
#Run verification process if desired
if verifyID:
    
    #Run ID with the .osim file
    
    #Set the path to generic templates
    pathGenericTemplates = os.path.join(baseDir, 'OpenSimPipeline')    
    pathGenericIDFolder = os.path.join(pathGenericTemplates, 'InverseDynamics')
    pathGenericIDSetupFile = os.path.join(pathGenericIDFolder, 'Setup_InverseDynamics.xml')
    
    #Create and run ID tool
    idTool = osim.InverseDynamicsTool(pathGenericIDSetupFile)
    idTool.setName('ID_withOsimAndIDTool')
    idTool.setModelFileName(outputModelFileName)
    idTool.setResultsDir(pathOutputExternalFunctionFolder)
    idTool.setCoordinatesFileName(os.path.join(pathGenericIDFolder, 'DefaultPosition_rajagopal.mot'))
    idTool.setOutputGenForceFileName('ID_withOsimAndIDTool.sto')
    pathSetupID = os.path.join(pathOutputExternalFunctionFolder, 'Setup_InverseDynamics.xml')
    idTool.printToXML(pathSetupID)
    idTool.run()
    
    #Extract torques from .osim + ID tool.    
    headers = []
    nCoordinatesAll = coordinateSet.getSize()
    for coord in range(nCoordinatesAll):                
        if (coordinateSet.get(coord).getName() == "pelvis_tx" or 
            coordinateSet.get(coord).getName() == "pelvis_ty" or 
            coordinateSet.get(coord).getName() == "pelvis_tz" or
            coordinateSet.get(coord).getName() == "knee_angle_r_beta" or 
            coordinateSet.get(coord).getName() == "knee_angle_l_beta"):
            suffix_header = "_force"
        else:
            suffix_header = "_moment"
        headers.append(coordinateSet.get(coord).getName() + suffix_header)
        
    ID_osim_df = storage_to_dataframe(os.path.join(
        pathOutputExternalFunctionFolder,"ID_withOsimAndIDTool.sto"), 
        headers)
    ID_osim = np.zeros((nCoordinates))
    count = 0
    for coordinate in coordinates:
        if (coordinate == "pelvis_tx" or 
            coordinate == "pelvis_ty" or 
            coordinate == "pelvis_tz"):
            suffix_header = "_force"
        else:
            suffix_header = "_moment"
        if 'beta' in coordinate:
            continue                
        ID_osim[count] = ID_osim_df.iloc[0][coordinate + suffix_header]
        count += 1
    
    #Extract torques from external function.
    import casadi as ca
    F_ext = '.dll'
    F = ca.external('F', os.path.join(
        pathOutputExternalFunctionFolder, externalFunctionName + F_ext))
    
    #Run checks for verification
    vec1 = np.zeros((nCoordinates*2, 1))
    vec1[::2, :] = 0.05   
    vec1[8, :] = -0.05
    vec2 = np.zeros((nCoordinates, 1))        
    if treadmill:
        vec4 = np.zeros((1, 1))
        vec3 = np.concatenate((vec1,vec2,vec4))
    else:            
        vec3 = np.concatenate((vec1,vec2))
    ID_F = (F(vec3)).full().flatten()[:nCoordinates]
    assert(np.max(np.abs(ID_osim - ID_F)) < 1e-6), (
        'error F vs ID tool {}'.format(np.max(np.abs(ID_osim - ID_F))))
    print('Verification torque generation: success')
    os.remove(os.path.join(
        pathOutputExternalFunctionFolder,"ID_withOsimAndIDTool.sto"))
    os.remove(os.path.join(
        pathOutputExternalFunctionFolder,"Setup_InverseDynamics.xml"))

# %% Input settings for simulation

#Note that the get_setup function comes from the UtilsDynamicSimulations folder 
#settings function, but won't import here for some reason, so it's defined earlier
# os.chdir(os.path.join(baseDir, 'UtilsDynamicSimulations'))
# sys.path.append(os.path.join(baseDir, 'UtilsDynamicSimulations'))
# from settingsOpenSimAD import get_setup
# os.chdir(pathMain)

#Get settings.
settings = get_setup(motion_type)

#Set the path to the motion file
pathMotionFile = os.path.join('data', 'Kinematics', trial_name + '.mot')

#Set time interval in settings        
settings['timeInterval'] = time_window

#Get demographics.    
settings['mass_kg'] = jaMass
settings['height_m'] = 1.75 ####made up - does this matter?

#Treadmill speed.
settings['treadmill_speed'] = treadmill_speed

#Trial name
settings['trial_name'] = trial_name

#OpenSim model name
settings['OpenSimModel'] = 'LaiUhlrich2022' #### is this needed?

# %% NOTE: UP TO HERE...EVERYTHING WORKING SO FAR...

"""

The above is the end of the processInputsOpenSimAD - so we now return to the
example_walking_opensimAD function where the settings are adjusted and constraints
added to the problem.

It may be useful at this point to step into an outside function and start using
the built in functions as they start to become more extensively used here...

"""























# %% ----- End of test_OpenSimAD_sprinting.py ----- %% #