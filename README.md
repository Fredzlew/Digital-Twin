
This README.txt file was generated on 30-07-2023 by Frederik Emil Serritzlew and Adam Alexander Hempler-Christiansen.


-------------------
GENERAL INFORMATION
-------------------

1. Author Information

     	1. Name: Frederik Emil Serritzlew
           Institution: Technical University of Denmark
           Address: Nils Koppels Allé, 2800 Kongens Lyngby, Denmark
           Email: s183309@dtu.dk

        2. Name: Adam Alexander Hempler-Christiansen
           Institution: Technical University of Denmark
           Address: Nils Koppels Allé, 2800 Kongens Lyngby, Denmark
           Email: s183535@dtu.dk


--------------------------
SHARING/ACCESS INFORMATION
-------------------------- 

1. Licenses/restrictions placed on the code: Non. 

---------------------
FOLDER & FILE OVERVIEW
---------------------

1. src - folder containing files for five sensor setup.
	1.1. data - folder coantining data files from experimental setup, simulated data and OMA data
		1.1.1. experimental_data - data files done on the experimental data
			1.1.1.1 Modal_par - data files coantining the modal parameters from OMA 
				1.1.1.1.1 FDDdamp_5_2_1.npy - damping ratio from EFDD for experimental data set with high damping
				1.1.1.1.2 FDDdamp_no_damp.npy - damping ratio from EFDD for experimental data set with no damping
				1.1.1.1.3 FDDfreq_5_2_1.npy - natural frequencies from EFDD for experimental data set with high damping
				1.1.1.1.4 FDDfreq_no_damp.npy - natural frequencies from EFDD for experimental data set with no damping
				1.1.1.1.5 FDDmodes_5_2_1.npy - mode shapes from EFDD for experimental data set with high damping
				1.1.1.1.6 FDDmodes_no_damp.npy - mode shapes from EFDD for experimental data set with no damping
				1.1.1.1.7 FDDPSD_5_2_1.npy - power spectral density from FDD for experimental data set with high damping
				1.1.1.1.8 FDDPSD_no_damp.npy - power spectral density from FDD for experimental data set with no damping
				1.1.1.1.9 FDDPSDfreq_5_2_1.npy - natural frequencies to the power spectral density from FDD for experimental data set with high damping
				1.1.1.1.10 FDDPSDfreq_no_damp.npy - natural frequencies to the power spectral density from FDD for experimental data set with no damping
				1.1.1.1.11 SSIdamp_5_2_1.npy - damping ratio from SSI-COV for experimental data set with high damping
				1.1.1.1.12 SSIdamp_no_damp.npy - damping ratio from SSI-COV for experimental data set with no damping
				1.1.1.1.13 SSIfreq_5_2_1.npy - natural frequencies from SSI-COV for experimental data set with high damping
				1.1.1.1.14 SSIfreq_no_damp.npy - natural frequencies from SSI-COV for experimental data set with no damping
				1.1.1.1.15 SSImodes_5_2_1.npy - mode shapes from SSI-COV for experimental data set with high damping
				1.1.1.1.16 SSImodes_no_damp.npy - mode shapes from SSI-COV for experimental data set with no damping
				1.1.1.1.17 SSIstab_5_2_1.npy - stabilization diagram from SSI-COV for experimental data set with high damping
				1.1.1.1.18 SSIstab_no_damp.npy - stabilization diagram from SSI-COV for experimental data set with no damping
			1.1.1.2 Data_de_dec_high.npy - detrending of the experimental data set with high damping 
			1.1.1.3 Data_de_dec_nodamp.npy - detrending of the experimental data set with no damping
			1.1.1.4 Filtdata_high.npy - Filtered of the experimental data set with high damping
			1.1.1.5 Filtdata_nodamp.npy - Filtered of the experimental data set with no damping
			1.1.1.6 Time_dec_high.npy - detrending the time of the experimental data set with high damping
			1.1.1.7 Time_dec_nodamp.npy - detrending the time of the experimental data set with no damping
			1.1.1.8 Time_high.npy - The time of the experimental data set with high damping
			1.1.1.9 Time_nodamp.npy - The time of the experimental data set with no damping
		1.1.2 simulated_data - folder able to coantain the simulation from the simalæted data and the modal parameters
			1.1.2.1 Modal_par - folder containing the modal parameters
				1.1.2.1.1 FDDdamp.npy - damping ratio from EFDD from simulated data set equal to 100 simulations
				1.1.2.1.2 FDDfreq.npy - natural frequencies from EFDD from simulated data set equal to 100 simulations
				1.1.2.1.3 FDDmodes.npy - mode shapes from EFDD from simulated data set equal to 100 simulations
				1.1.2.1.4 FDDomegas.npy - natural angular frequencies from EFDD from simulated data set equal to 100 simulations
				1.1.2.1.5 FDDPSD.npy - power sepctral density from FDD from 1 simulated data set 
				1.1.2.1.6 FDDPSDfreq.npy - natural frequencies to power sepctral density from FDD from simulated data set equal to 100 simulations
				1.1.2.1.7 SSI_stab_sim.npy - stabilization diagram from SSI-COV from 1 simulated data set 
				1.1.2.1.8 SSIcovDamp.npy - damping ratio from SSSI-COV from simulated data set equal to 100 simulations
				1.1.2.1.9 SSIcovfreq.npy - natural frequencies from SSI-COV from simulated data set equal to 100 simulations
				1.1.2.1.1 SSIcovModes.npy - mode shapes from SSI-COV from simulated data set equal to 100 simulations
				1.1.2.1.1 SSIcovOmegas.npy - natural angular frequencies from SSI-COV from simulated data set equal to 100 simulations
				1.1.2.1.1 SSIdatdamp.npy - damping ratio from SSI-DAT from simulated data set equal to 100 simulations
				1.1.2.1.1 SSIdatfreq.npy - natural frequencies from SSI-DAT from simulated data set equal to 100 simulations
				1.1.2.1.1 SSIdatmodes.npy - mode shapes from SSI-DAT from simulated data set equal to 100 simulations
		1.1.3 modelprop.mat - model properties and modal parameters from the FE model
	1.2 Fatigue - folder containing files used for fatigue evaluation
		1.2.1 Rainflow_data - folder containing results from rainflow counting
			1.2.1.1 me_rainflow_high.csv - results from rainflow counting for measured data with high damping
			1.2.1.2 me_rainflow_nodamp.csv - results from rainflow counting for measured data with no damping	
			1.2.1.3 pr_rainflow_high.csv - results from rainflow counting for predicted data with high damping
			1.2.1.4 pr_rainflow_nodamp.csv - results from rainflow counting for predicted data with no damping
		1.2.2 Rainflow.py - file that performs rainflow counting
	1.3 Model_updating - folder containing files used for model updating
		1.3.1 Costfunctions - folder with files used for model updating based on cost functions
			1.3.1.1 data_updated_par - folder containing updated modal parameters after model updating
				1.3.1.1.1 SSIfreq_high.mat - modal parameters after eigenvalue cost function with high damping
				1.3.1.1.2 SSIfreq_no_damp.mat - modal parameters after eigenvalue cost function with no damping
				1.3.1.1.3 SSIfreqmodeEIL_high.mat - modal parameters after eigenvalue+modeshape cost function with high damping and geometrical/material properties as optimization parameters
				1.3.1.1.4 SSIfreqmodeEIL_no_damp.mat - modal parameters after eigenvalue+modeshape cost function with no damping and geometrical/material properties as optimization parameters
				1.3.1.1.5 SSIfreqmode_high.mat - modal parameters after eigenvalue+modeshape cost function with high damping
				1.3.1.1.6 SSIfreqmode_no_damp.mat - modal parameters after eigenvalue+modeshape cost function with no damping
				1.3.1.1.7 SSImode_high.mat - modal parameters after modeshape cost function with high damping
				1.3.1.1.8 SSImode_mac_high.mat - modal parameters after modeshape with MAC cost function with high damping
				1.3.1.1.9 SSImode_mac_no_damp.mat - modal parameters after modeshape with MAC cost function with no damping
				1.3.1.1.10 SSImode_no_damp.mat - modal parameters after modeshape cost function with no damping
			1.3.1.2 functions - folder containing cost functions and other relevant functions
				1.3.1.2.1 costfunSSIfreq.m - function that contains the eigenvalue cost function
				1.3.1.2.2 costfunSSIfreqmode.m - function that contains the eigenvalue+modeshape cost function
				1.3.1.2.3 costfunSSIfreqmodeEIL.m - function that contains the eigenvalue+modeshape cost function with geometrical/material properties as optimization parameters
				1.3.1.2.4 costfunSSImode.m - function that contains the modeshape cost function
				1.3.1.2.5 costfunSSImode_mac.m - function that contains the modeshape with MAC cost function
				1.3.1.2.6 crossMAC.m - function to calculate the CrossMAC value of two sets of vectors
				1.3.1.2.7 getGlobalx.m - function that gets a global variable
				1.3.1.2.8 setGlobalx.m - function that defines a global variable
			1.3.1.3 FE_model_update.m - file that performs the model updating 
			1.3.1.4 Plt_cost_high.m - file that plots the results with high damping
			1.3.1.5 Plt_cost_no_damp.m - file that plots the results with no damping
		1.3.2 Sensitivity_method - folder with files used for model updating based on the sensitivity method
			1.3.2.1 data_updated_par_sens - folder containing results from model updating peformed with the sensitivity method
				1.3.2.1.1 Eigenvalue_Mode_shape_residual_L_curve_high.mat - file containing results from L-curve of data with high damping based on eigenvalue+modeshape residual
				1.3.2.1.2 Eigenvalue_Mode_shape_residual_L_curve_no_damp.mat - file containing results from L-curve of data with no damping based on eigenvalue+modeshape residual
				1.3.2.1.3 Eigenvalue_Mode_shape_residual_high.mat - file containing the updated modal parameters and stiffnesses for the data with high damping based on eigenvalue+modeshape residual
				1.3.2.1.4 Eigenvalue_Mode_shape_residual_no_damp.mat - file containing the updated modal parameters and stiffnesses for the data with no damping based on eigenvalue+modeshape residual
				1.3.2.1.5 Eigenvalue_residual_L_curve_high.mat - file containing results from L-curve of data with high damping based on eigenvalue residual
				1.3.2.1.6 Eigenvalue_residual_L_curve_no_damp.mat - file containing results from L-curve of data with no damping based on eigenvalue residual
				1.3.2.1.7 Eigenvalue_residual_high.mat - file containing the updated modal parameters and stiffnesses for the data with high damping based on eigenvalue residual
				1.3.2.1.8 Eigenvalue_residual_no_damp.mat - file containing the updated modal parameters and stiffnesses for the data with no damping based on eigenvalue residual
				1.3.2.1.9 Mode_shape_residual_L_curve_high.mat - file containing results from L-curve of data with high damping based on modeshape residual
				1.3.2.1.10 Mode_shape_residual_L_curve_no_damp.mat - file containing results from L-curve of data with no damping based on modeshape residual
				1.3.2.1.11 Mode_shape_residual_high.mat - file containing the updated modal parameters and stiffnesses for the data with high damping based on modeshape residual
				1.3.2.1.12 Mode_shape_residual_no_damp.mat - file containing the updated modal parameters and stiffnesses for the data with no damping based on modeshape residual
			1.3.2.2 mot_ex - folder containing the calculations for the examples given in "The sensitivity method in finite element model updating: A tutorial" by John E. Mottershead et al. 2011.
				1.3.2.2.1 Case1_eigen_residual.m - file containing an example of implementation of the eigenvalue residual for the examples
				1.3.2.2.2 Cases.m - file containing the calculations of the examples
			1.3.2.3 Eigenvalue_Mode_shape_residual.m - file that performs model updating with the eigenvalue+modeshape residual
			1.3.2.4 Eigenvalue_Mode_shape_residual_stffmass.m - test-file that performs model updating with the eigenvalue+modeshape residual with both stiffnesses and masses
			1.3.2.5 Eigenvalue_residual.m - file that performs model updating with the eigenvalue residual
			1.3.2.6 Mode_shape_residual.m - file that performs model updating with the modeshape residual
			1.3.2.7 Plt_sens_method_high.m - file that plots the results with high damping
			1.3.2.8 Plt_sens_method_no_damp.m - file that plots the results with no damping
		1.3.3 Plt_compare_method_high.m - file that plots all the results from model updating with high damping
		1.3.4 Plt_compare_method_no_damp.m - file that plots all the results from model updating with no damping
	1.4 npy-matlab-master - folder containing scripts enabling python files to read in matlab. Since we did not create this folder ourselves, this will not be explained any further. A specific ReadMe file can be found in the folder
	1.5 OMA - folder containing the different OMA scripts
		1.5.1 Experimental - the OMA scripts used on experimental data
			1.5.1.1 FDD_exp.ipy - OMA script that runs FDD and EFDD on the experimental data 
			1.5.1.2 SSI_modeshapes.m - script that show the mode shapes from the SSI-COV modal parameters
			1.5.1.3	SSIcov_exp.ipy - OMA script that runs SSI-COV on the experimental data 
		1.5.2 Simulated - the OMA scripts used on the simulated data
			1.5.2.1 FDD_sim.ipy - OMA script that runs FDD and EFDD on the simulated data
			1.5.2.2 FDD_sim_PSD.ipy - OMA script that runs FDD to obtain the power spectral density plot on the simulated data
			1.5.2.3 Simulated_data_compare.m - script that compare all the 3 OMA methods
			1.5.2.4 SSIcov_sim.ipy -  OMA script that runs SSI-COV on the simulated data
			1.5.2.5 SSIcov_sim_stab.ipy - OMA script that runs SSI-COV to obtain the stabilization diagram on the simulated data
			1.5.2.6 SSIdat_sim.ipy -  OMA script that runs SSI-DAT on the simulated data
		1.5.3 PSD_Plot_sim.m - script that plots the power spectral density for the simulated data
		1.5.4 Stabilization_diagram.m - script that plots the stabilization diagram for experimental data
		1.5.5 Stabilitzation_diagram_sim.m - script that plots the stabilization diagram for the simulated data
	1.6 Virtual_sensing - folder that contain the scripts that run virtual sensing and detrending and filtering of the data
		1.6.1 Experimental - folder that contains the scripts done on the experimental data		
			1.6.1.1 Filtered_data -	folder coanting the data gain from 1.6.1.2 Multiband.py 
				1.6.1.1.1 data_filt_all_no_damp.npy - experimental data set with no damping 
				1.6.1.1.2 data_filt_first_modes_no_damp.npy - file containing a low-pass filter on the experimental data set with no damping 
				1.6.1.1.3 data_filt_mode3_no_damp.npy - file containing a band-pass filter on the experimental data set with no damping for the 3 mode
				1.6.1.1.4 data_filt_mode4_no_damp.npy - file containing a band-pass filter on the experimental data set with no damping for the 4 mode
				1.6.1.1.5 data_filt_mode5_no_damp.npy - file containing a band-pass filter on the experimental data set with no damping for the 5 mode
				1.6.1.1.6 data_filt_mode1_no_damp.npy - file containing a low-pass filter on the experimental data set with no damping for the 1 mode
				1.6.1.1.7 data_filt_mode2_no_damp.npy - file containing a band-pass filter on the experimental data set with no damping for the 2 mode
				1.6.1.1.8 data_filt_0_cut_no_damp.npy - file containing the 2 first frequencies bands for the experimental data set with no damping 
				1.6.1.1.9 data_filt_cut_end_no_damp.npy - file containing the 3 last frequencies bands for the experimental data set with no damping 
				1.6.1.1.10 data_filt_no_damp.npy - experimental data set with no damping that has a cut off the last mode
				1.6.1.1.11 data_predicted_nodamp.txt - file containing the predicted displacements from the multiband with 5 bands with the experimatal data set with high damping
				1.6.1.1.12 data_predicted_high.txt - file containing the predicted displacements from the multiband with 5 bands with the experimatal data set with no damping
				1.6.1.1.13 data_filt_all_high.npy - experimental data set with high damping 
				1.6.1.1.14 data_filt_first_modes_high.npy - file containing a low-pass filter on the experimental data set with high damping 
				1.6.1.1.15 data_filt_mode3_high.npy - file containing a band-pass filter on the experimental data set with high damping for the 3 mode
				1.6.1.1.16 data_filt_mode4_high.npy - file containing a band-pass filter on the experimental data set with high damping for the 4 mode
				1.6.1.1.17 data_filt_mode5_high.npy - file containing a band-pass filter on the experimental data set with high damping for the 5 mode
				1.6.1.1.18 data_filt_mode1_high.npy - file containing a band-pass filter on the experimental data set with high damping for the 1 mode
				1.6.1.1.19 data_filt_mode2_high.npy - file containing a band-pass filter on the experimental data set with high damping for the 2 mode
				1.6.1.1.20 data_filt_0_cut_high.npy - file containing the 2 first frequencies bands for the experimental data set with high damping 
				1.6.1.1.21 data_filt_cut_end_high.npy - file containing the 3 last frequencies bands for the experimental data set with high damping 
				1.6.1.1.22 data_filt_high.npy - experimental data set with high damping that has a cut off the last mode
			1.6.1.2 Multiband.py - script that perfom filtering on the experimental data and divided it up in different bands	
			1.6.1.3 Virtual_sensing_all_data.m - script that runs the virtual sensing on experimental data which is not filtered
			1.6.1.4 Virtual_sensing_first_modes.m -	script that runs the virtual sensing on experimental data which is filtered with a low-pass filter
			1.6.1.5 Virtual_sensing_multi.m - script that runs the virtual sensing on experimental data which is filtered with a low pass and a band-pass	
			1.6.1.6 Virtual_sensing_multi_5_modes.m - script that runs the virtual sensing with 5 bands on experimental data which is filtered with a low-pass filter and four band-passes	
			1.6.1.7 Virtual_sensing_multi_5_modes_different_mode_shapes.m -	script that runs virtual sensing with 5 bands on experimental data with all the data from the different modal updating method 
		1.6.2 function - folder containing files used for virtual sensing
            		1.6.2.1 VirtualSensVal.m - function that peforms virtual sensing based on the MDE method
        	1.6.3 Simulated - folder containing files used for virtual sensing of the simulated data
            		1.6.3.1 Filtered_data_sim - folder containing the filtered data
                		1.6.3.1.1 data_filt_0_cut_sim.npy - file containing the filtered data in the first frequency band
                		1.6.3.1.2 data_filt_cut_end_sim.npy - file containing the filtered data in the second frequency band
                		1.6.3.1.3 data_filtdata_all_sim.npy - file containing the filtered data 
                		1.6.3.1.4 data_filtdata_first_modes_sim.npy - file containing the filtered data where higher modes are excluded
                		1.6.3.1.5 sim_data.npy - file coanting the simulated data after detrending
            		1.6.3.2 simmultiband.py - file splitting the data into two frequency bands
            		1.6.3.3 Virtual_sensing_all_data_sim.m - file used to perform virtual sensing on all the data
            		1.6.3.4 Virtual_sensing_first_modes_sim.m - file used to perform virtual sensing where the higher modes are excluded
            		1.6.3.5 Virtual_sensing_multi_sim.m - file used to perform virtual sensing with two frequency bands
    1.7 Data_simulation.m - file used to simulate data
    1.8 Exp_data.ipy - file comparing raw experimental data to the data after filtering, detrending and down-sampling has been applied
    1.9 FE_model.m - file containing the FE-model

2 src_operational - folder containing files for three sensor setup.
	2.1 data - folder coantining data files from experimental setup, simulated data and OMA data
		2.1.1 experimental_data - data files done on the experimental data
			2.1.1.1 Modal_par_3_sensors - data files coantining the modal parameters from OMA 
				2.1.1.1.1 SSIdamp_5_2_1.npy - damping ratio from SSI-COV for experimental data set with high damping
				2.1.1.1.2 SSIdamp_no_damp.npy - damping ratio from SSI-COV for experimental data set with no damping
				2.1.1.1.3 SSIfreq_5_2_1.npy - natural frequencies from SSI-COV for experimental data set with high damping
				2.1.1.1.4 SSIfreq_no_damp.npy - natural frequencies from SSI-COV for experimental data set with no damping
				2.1.1.1.5 SSImodes_5_2_1.npy - mode shapes from SSI-COV for experimental data set with high damping
				2.1.1.1.6 SSImodes_no_damp.npy - mode shapes from SSI-COV for experimental data set with no damping
				2.1.1.1.7 SSIstab_5_2_1.npy - stabilization diagram from SSI-COV for experimental data set with high damping
				2.1.1.1.8 SSIstab_no_damp.npy - stabilization diagram from SSI-COV for experimental data set with no damping
		2.1.3 modelprop.mat - model properties and modal parameters from the FE model
	2.2 Fatigue - folder containing files used for fatigue evaluation
		2.2.1 Rainflow_data_3_sensors - folder containing results from rainflow counting
			2.2.1.1 me_rainflow_3_sensors_high.csv - results from rainflow counting for measured data with high damping
			2.2.1.2 me_rainflow_3_sensors_nodamp.csv - results from rainflow counting for measured data with no damping	
			2.2.1.3 pr_rainflow_3_sensors_high.csv - results from rainflow counting for predicted data with high damping
			2.2.1.4 pr_rainflow_3_sensors_nodamp.csv - results from rainflow counting for predicted data with no damping
		2.2.2 Rainflow_3_sensors.py - file that performs rainflow counting
	2.3 Model_updating - folder containing files used for model updating
		2.3.1 Costfunctions - folder with files used for model updating based on cost functions
			2.3.1.1 data_updated_par_3_sensors - folder containing updated modal parameters after model updating
				2.3.1.1.1 SSIfreq_3_sensors_high.mat - modal parameters after eigenvalue cost function with high damping
				2.3.1.1.2 SSIfreq_3_sensors_no_damp.mat - modal parameters after eigenvalue cost function with no damping
				2.3.1.1.3 SSIfreqmodeEIL_3_sensors_high.mat - modal parameters after eigenvalue+modeshape cost function with high damping and geometrical/material properties as optimization parameters
				2.3.1.1.4 SSIfreqmodeEIL_3_sensors_no_damp.mat - modal parameters after eigenvalue+modeshape cost function with no damping and geometrical/material properties as optimization parameters
				2.3.1.1.5 SSIfreqmode_3_sensors_high.mat - modal parameters after eigenvalue+modeshape cost function with high damping
				2.3.1.1.6 SSIfreqmode_3_sensors_no_damp.mat - modal parameters after eigenvalue+modeshape cost function with no damping
				2.3.1.1.7 SSImode_3_sensors_high.mat - modal parameters after modeshape cost function with high damping
				2.3.1.1.8 SSImode_mac_3_sensors_high.mat - modal parameters after modeshape with MAC cost function with high damping
				2.3.1.1.9 SSImode_mac_3_sensors_no_damp.mat - modal parameters after modeshape with MAC cost function with no damping
				2.3.1.1.10 SSImode_3_sensors_no_damp.mat - modal parameters after modeshape cost function with no damping
			2.3.1.2 functions - folder containing cost functions and other relevant functions
				2.3.1.2.1 costfunSSIfreq_3_sensors.m - function that contains the eigenvalue cost function
				2.3.1.2.2 costfunSSIfreqmode_3_sensors.m - function that contains the eigenvalue+modeshape cost function
				2.3.1.2.3 costfunSSIfreqmodeEIL_3_sensors.m - function that contains the eigenvalue+modeshape cost function with geometrical/material properties as optimization parameters
				2.3.1.2.4 costfunSSImode_3_sensors.m - function that contains the modeshape cost function
				2.3.1.2.5 costfunSSImode_ma_3_sensorsc.m - function that contains the modeshape with MAC cost function
				2.3.1.2.6 crossMAC.m - function to calculate the CrossMAC value of two sets of vectors
				2.3.1.2.7 getGlobalx.m - function that gets a global variable
				2.3.1.2.8 setGlobalx.m - function that defines a global variable
			2.3.1.3 FE_model_3_sensors_update.m - file that performs the model updating 
			2.3.1.4 Plt_cost_3_sensors_high.m - file that plots the results with high damping
			2.3.1.5 Plt_cost_3_sensors_no_damp.m - file that plots the results with no damping
		2.3.2 Sensitivity_method - folder with files used for model updating based on the sensitivity method
			2.3.2.1 data_updated_par_sens_3_sensors - folder containing results from model updating peformed with the sensitivity method
				2.3.2.1.1 Eigenvalue_Mode_shape_residual_L_curve_3_sensors_high.mat - file containing results from L-curve of data with high damping based on eigenvalue+modeshape residual
				2.3.2.1.2 Eigenvalue_Mode_shape_residual_L_curve_3_sensors_no_damp.mat - file containing results from L-curve of data with no damping based on eigenvalue+modeshape residual
				2.3.2.1.3 Eigenvalue_Mode_shape_residual_3_sensors_high.mat - file containing the updated modal parameters and stiffnesses for the data with high damping based on eigenvalue+modeshape residual
				2.3.2.1.4 Eigenvalue_Mode_shape_residual_3_sensors_no_damp.mat - file containing the updated modal parameters and stiffnesses for the data with no damping based on eigenvalue+modeshape residual
				2.3.2.1.5 Eigenvalue_residual_L_curve_3_sensors_high.mat - file containing results from L-curve of data with high damping based on eigenvalue residual
				2.3.2.1.6 Eigenvalue_residual_L_curve_3_sensors_no_damp.mat - file containing results from L-curve of data with no damping based on eigenvalue residual
				2.3.2.1.7 Eigenvalue_residual_3_sensors_high.mat - file containing the updated modal parameters and stiffnesses for the data with high damping based on eigenvalue residual
				2.3.2.1.8 Eigenvalue_residual_3_sensors_no_damp.mat - file containing the updated modal parameters and stiffnesses for the data with no damping based on eigenvalue residual
				2.3.2.1.9 Mode_shape_residual_L_curve_3_sensors_high.mat - file containing results from L-curve of data with high damping based on modeshape residual
				2.3.2.1.10 Mode_shape_residual_L_curve_3_sensors_no_damp.mat - file containing results from L-curve of data with no damping based on modeshape residual
				2.3.2.1.11 Mode_shape_residual_3_sensors_high.mat - file containing the updated modal parameters and stiffnesses for the data with high damping based on modeshape residual
				2.3.2.1.12 Mode_shape_residual_3_sensors_no_damp.mat - file containing the updated modal parameters and stiffnesses for the data with no damping based on modeshape residual
			2.3.2.3 Eigenvalue_Mode_shape_residual_3_sensors.m - file that performs model updating with the eigenvalue+modeshape residual
			2.3.2.4 Eigenvalue_residual_3_sensors.m - file that performs model updating with the eigenvalue residual
			2.3.2.5 Mode_shape_residual_3_sensors.m - file that performs model updating with the modeshape residual
			2.3.2.6 Plt_sens_method_3_sensors_high.m - file that plots the results with high damping
			2.3.2.7 Plt_sens_method_3_sensors_no_damp.m - file that plots the results with no damping
		2.3.3 Plt_compare_method_high.m - file that plots all the results from model updating with high damping
		2.3.4 Plt_compare_method_no_damp.m - file that plots all the results from model updating with no damping
	2.4 npy-matlab-master - folder containing scripts enabling python files to read in matlab. Since we did not create this folder ourselves, this will not be explained any further. A specific ReadMe file can be found in the folder
	2.5 OMA - folder containing the different OMA scripts
		2.5.1 Experimental - the OMA scripts used on experimental data
			2.5.1.1 FDD_exp_3_sensors.ipy - OMA script that runs FDD and EFDD on the experimental data 
			2.5.1.2 SSI_modeshapes.m - script that show the mode shapes from the SSI-COV modal parameters
			2.5.1.3	SSIcov_exp_3_sensors.ipy - OMA script that runs SSI-COV on the experimental data 
	2.6 Virtual_sensing - folder that contain the scripts that run virtual sensing and detrending and filtering of the data
		2.6.1 Experimental - folder that contains the scripts done on the experimental data		
			2.6.1.1 Filtered_data -	folder coanting the data gain from 1.6.1.2 Multiband.py 
				2.6.1.1.1 data_filt_all_3_sensors_no_damp.npy - experimental data set with no damping 
				2.6.1.1.2 data_filt_first_modes_3_sensors_no_damp.npy - file containing a low-pass filter on the experimental data set with no damping 
				2.6.1.1.3 data_filt_mode3_3_sensors_no_damp.npy - file containing a band-pass filter on the experimental data set with no damping for the 3 mode
				2.6.1.1.4 data_filt_mode4_3_sensors_no_damp.npy - file containing a band-pass filter on the experimental data set with no damping for the 4 mode
				2.6.1.1.5 data_filt_mode5_3_sensors_no_damp.npy - file containing a band-pass filter on the experimental data set with no damping for the 5 mode
				2.6.1.1.6 data_filt_mode1_3_sensors_no_damp.npy - file containing a low-pass filter on the experimental data set with no damping for the 1 mode
				2.6.1.1.7 data_filt_mode2_3_sensors_no_damp.npy - file containing a band-pass filter on the experimental data set with no damping for the 2 mode
				2.6.1.1.8 data_filt_0_cut_3_sensors_no_damp.npy - file containing the 2 first frequencies bands for the experimental data set with no damping 
				2.6.1.1.9 data_filt_cut_end_3_sensors_no_damp.npy - file containing the 3 last frequencies bands for the experimental data set with no damping 
				2.6.1.1.10 data_filt_3_sensors_no_damp.npy - experimental data set with no damping that has a cut off the last mode
				2.6.1.1.11 data_predicted_3_sensors_nodamp.txt - file containing the predicted displacements from the multiband with 5 bands with the experimatal data set with high damping
				2.6.1.1.12 data_predicted_3_sensors_high.txt - file containing the predicted displacements from the multiband with 5 bands with the experimatal data set with no damping
				2.6.1.1.13 data_filt_all_3_sensors_high.npy - experimental data set with high damping 
				2.6.1.1.14 data_filt_first_modes_3_sensors_high.npy - file containing a low-pass filter on the experimental data set with high damping 
				2.6.1.1.15 data_filt_mode3_3_sensors_high.npy - file containing a band-pass filter on the experimental data set with high damping for the 3 mode
				2.6.1.1.16 data_filt_mode4_3_sensors_high.npy - file containing a band-pass filter on the experimental data set with high damping for the 4 mode
				2.6.1.1.17 data_filt_mode5_3_sensors_high.npy - file containing a band-pass filter on the experimental data set with high damping for the 5 mode
				2.6.1.1.18 data_filt_mode1_3_sensors_high.npy - file containing a band-pass filter on the experimental data set with high damping for the 1 mode
				2.6.1.1.19 data_filt_mode2_3_sensors_high.npy - file containing a band-pass filter on the experimental data set with high damping for the 2 mode
				2.6.1.1.20 data_filt_0_cut_3_sensors_high.npy - file containing the 2 first frequencies bands for the experimental data set with high damping 
				2.6.1.1.21 data_filt_cut_end_3_sensors_high.npy - file containing the 3 last frequencies bands for the experimental data set with high damping 
				2.6.1.1.22 data_filt_3_sensors_high.npy - experimental data set with high damping that has a cut off the last mode
			2.6.1.2 Multiband_3_sensors.py - script that perfom filtering on the experimental data and divided it up in different bands	
			2.6.1.3 Virtual_sensing_all_data_3_sensors.m - script that runs the virtual sensing on experimental data which is not filtered
			2.6.1.4 Virtual_sensing_first_modes_3_sensors.m - script that runs the virtual sensing on experimental data which is filtered with a low-pass filter
			2.6.1.5 Virtual_sensing_multi_3_sensors.m - script that runs the virtual sensing on experimental data which is filtered with a low pass and a band-pass	
			2.6.1.6 Virtual_sensing_multi_5_modes_3_sensors.m - script that runs the virtual sensing with 5 bands on experimental data which is filtered with a low-pass filter and four band-passes	
		2.6.2 function - folder containing files used for virtual sensing
            		2.6.2.1 VirtualSensVal_3_sensors.m - function that peforms virtual sensing based on the MDE method
    2.7 FE_model_3_sensors.m - file containing the FE-model


-----------------
EXPERIMENTAL DATA
-----------------

The experimental data is available at https://data.mendeley.com/datasets/77y99kvt9k/1.
For this thesis, a data set with no damping "data.txt" and a data set with high damping 
"data_5_2_1.txt" have been used. The experimental data sets should be saved at:
*src->data->experimental_data (1.1.1).

For more information on the experimental data, the user is advised to:
Bajric, A., Høgsberg, J., Identification of damping and complex modes in structural vibrations, 2017.

--------------------------
METHODOLOGICAL INFORMATION
--------------------------
The following will be an overview of how to run the files. This will be valid for both
the setup with five and three sensors. Five sensors are used in "src", while three sensors
are used in "src_operational". The chronology will be presented for the five sensor setup,
with experimental data.

Note: Different OMA techniques, cost functions or virtual sensing schemes are available. 
Simply choose a different folder or file.
Note: By running the simulated python scripts it is important to import your own path in the start of the python scripts.

1. Run "FE_model" to obtain the finite element model of the structure.
*src->FE_model.m (1.9)

2. Run "Data_simulation" to simulate data (not necessary, if only experimental data is considered).
*src->Data_simulation.m (1.7)

3. Run "SSIcov_exp" to perform OMA.
*src->OMA->Experimental->SSIcov_exp.ipy (1.5.1.1)

4. Run "Eigenvalue_Mode_shape_residual" to perform model updating with the sensitivity method.
*src->Model_updating->Sensitivity_method->Eigenvalue_Mode_shape_residual.m (1.3.2.3)

5. Run "Multiband" to divide the data into frequency bands.
*src->Virtual_sensing->Experimental->Multiband.py (1.6.1.2)

6. Run "Virtual_sensing_multi_5_modes" to perform virtual sensing with five frequency bands.
*src->Virtual_sensing->Experimental->Virtual_sensing_multi_5_modes.m (1.6.1.6)

7. Run "Rainflow" to perform rainflow counting on measured and predicted data.
*src->Fatigue->Rainflow.py (1.2.2)

