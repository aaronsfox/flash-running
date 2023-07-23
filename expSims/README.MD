This folder contains a subset of the experimental data provided by Dorn et al. (2012). Specifically, the scaled OpenSim model (see note below), the static calibration trial ('JA1Static05.c3d'), and the and 9.0 metres per second sprinting trial ('JA1Gait35.c3d' — relabelled to 'JA1Gait35_9ms.c3d').

The original model has been slightly updated to allow processing, likely due to some inconsistencies with OpenSim 4.0. The muscles in the force set were removed as these had some parameters that didn't agree with OpenSim, and the model was updated to use 4.0 coordinate notation ('JA1_SCALED_Osim40.osim').

The Python scripts included in an adjacent **code** folder convert this data to an OpenSim format ('convertData.py') and run a tracking simulation on the data ('processData.py') to generate dynamically consistent sprinting data as a guide for further simulations.

## References

Dorn et al. (2012). Muscular strategy shift in human running: dependence of running speed on hip and ankle muscle performance. *J Exp Biol*, 215: 1944-1956.