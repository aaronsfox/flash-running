This folder contains a subset of the experimental data provided by Dorn et al. (2012). Some modifications have been made and certain files renamed for simplicity. A summary of the files is as follows:

- **processRealWord.py:** This is the main data processing script within this folder. It contains all the necessary code to process the data from inverse kinematics through to the simulations. See the comments within this script for more details on the code.
- **data/convertData.py:** This is a supplementary Python script that converts the original Dorn et al. (2012) data files to what is outlined below with respect to the 'data' folder. There are certain dependencies required to run this script (e.g. OpenSim, btk) - where these can be found by examining the packages imported in the script.
- **models/base3DModel.osim:** This was originally labelled 'JA1_SCALED.osim' within the Dorn et al. (2012) dataset. Some modifications to the model were made to ensure compatibility with OpenSim 4.3, namely the muscles were removed and the model was updated to use OpenSim 4.x joint system notation. A basic set of contact spheres were added to the bottom of the feet to be used in the later tracking & predictive simulations.
- **models/gait9dof18musc_Ong_et_al_Moco.osim:** This is a 2D model with a simplified muscle set from Ong et al. (2019). This model isn't used in any simulations, but the lower limb muscle set from this model will be applied to a 2D model created from the base3DModel to be used in simulations.
- **data/static.trc:** This is the marker data from the static trial (i.e. 'JA1Static05.c3d') included with the Dorn et al. (2012) dataset. The marker data from the .c3d file were extracted, rotated to the OpenSim coordinate system, and converted to the .trc format.
- **data/sprint.trc:** This is the marker data from the sprint trial (i.e. 'JA1Gait35.c3d') included with the Dorn et al. (2012) dataset. The marker data from the .c3d file were extracted, rotated to the OpenSim coordinate system, and converted to the .trc format.
- **data/sprint_grf.mot:** This is the GRF data from the sprint trial (i.e. 'JA1Gait35.c3d') included with the Dorn et al. (2012) dataset. The force data from the .c3d file were extracted, rotated to the OpenSim coordinate system, and converted to the .mot format.
- **data/sprint_grf.xml:** This is the external loads .xml file created to align with the 'sprint_grf.mot' file.



## References

Dorn et al. (2012). Muscular strategy shift in human running: dependence of running speed on hip and
ankle muscle performance. *J Exp Biol*, 215: 1944-1956.

Ong et al. (2019). Predicting gait adaptations due to ankle plantarflexor muscle weakness and contracture using physics-based musculoskeletal simulations. *PLoS Comput Biol*, 15: e1006993.