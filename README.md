# csi_project_AGV
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Realization and simulation of an AGV dynamic model for the design of H-infinite and LQG optimal controllers.
The project has been realised for the course of Uncertain Systems Control at University of Pisa (Course Of Studies in Robotics).

# Requirements
MATLAB 2022b or later versions.

# Instructions
Inside the FINAL_PROJECT folder:

Execute the script "LQG_no_actuators.m" to design a LQG optimal controller without actuator transfer functions. For the simulation execute the simulink model "LQG_no_integrator.slx" for linearized system and the simulink model "LQG_nonlinear.slx" for non linear system.

Execute the script "LQG_integrator.m" to design a LQG optimal controller without actuator transfer functions and with integral action.

Execute the script "LQG_act.m" to design a LQG optimal controller with actuator transfer functions. For the simulation execute the simulink model "LQG_no_integrator_act.slx" for linearized system.

Execute the script "Hinf_controller_hinfsyn.m" and "Hinf_controller_mixsyn.m" to design a H-infinite controller with two differente methods. For the simulation execute the simulink model "mixed_sensitivity.slx" for linearized system and the simulink model "mixed_sensitivity_nonlin.m" for non linear system.

Execute the script "DK_iteration.m" to design an H-infinite controller with DK-iteration. For the simulation execute the simulink model "dk_lin.slx" for linearized system and the simulink model "dk_non_lin.m" for non linear system.

Inside the LFTDATA folder:

Execute the script "Hinf_controller_Hinf.m" and "Hinf_controller_mixsyn.m" to design a H-infinite controller with two differente methods. For the simulation execute the simulink model "mixed_sensitivityhinf.slx" for linearized system.

Execute the script "DK_iteration.m" to design an H-infinite controller with DK-iteration. For the simulation execute the simulink model "dk.slx" for linearized system and the simulink model "dk_nonlin.m" for non linear system.
