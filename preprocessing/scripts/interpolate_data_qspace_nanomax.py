
# The script converts the 5D scan aligned-data to orthogonal q-space for evaluating the Bragg peaks visually
# 
import os
import numpy as np
import sys
sys.path.append('/home/dzhigd/Software/') # Important to check!
# import nanomax_tools.preprocessing.qspace_utils as qu

# INPUTS ######################################################################
year        = "2020"                                                           #The year for the experiemnt
beamtimeID  = "2020101408"                                                     #The beamtimeID
sample_name = r"0002_sample_P246_AB"                                           #The name for the p10 newfile
scan_number = np.linspace(109,128,20)                                          #The scan numbers, can be the list
rocking_motor = "gontheta" # "gonphi"

# Processing path
processing_root = '/home/dzhigd/work/projects/Qdevs_2020_NanoMAX/data/scan_%d_%d'%(scan_number[0],scan_number[-1])

# Data path
path_root = '/media/dzhigd/My Passport/DDzhigaev/Data/MAXIV/NanoMax/%s/process/%s_%d_%d/data_aligned'%(beamtimeID,sample_name,scan_number[0],scan_number[-1])

for ii in range(0,len(scan_number)):
	print(('\033[1;33;40m #### Processing scan %d')%(int(scan_number[ii])))
	data_npz_path = os.path.join(path_root, "scan_%06d_merlin.npz"%(int(scan_number[ii])))
    # Load the data 
	input_data = np.load(data_npz_path,allow_pickle=True)
	dictionary = input_data['arr_0'].item()

	data_xrd = dictionary['data_xrd']
	xrd_roi  = dictionary['xrd_roi']
	scan_position_x = dictionary['scan_position_x']
	scan_position_y = dictionary['scan_position_y']
	scan_position_z = dictionary['scan_position_z']
	motor_positions = dictionary['motor_positions']
	rocking_motor   = dictionary['rocking_motor']
	rocking_angles  = dictionary['rocking_angles']

	data_t = np.sum(np.sum(data_xrd,axis=0),axis=0)

	if ii == 0:
		data = np.empty((len(scan_number), xrd_roi[1]-xrd_roi[0],xrd_roi[3]-xrd_roi[2]))

		if rocking_motor == "gontheta":		
			gontheta = np.empty((len(scan_number)))
		elif rocking_motor == "gonphi":
			gonphi = np.empty((len(scan_number)))			
		else:
			print('-- Unknown rocking motor! Saved Nothing!')
	
	data[ii,:,:] = data_t

	if rocking_motor == "gontheta":		
		gontheta[ii] = rocking_angles
	elif rocking_motor == "gonphi":		
		gonphi[ii] = rocking_angles

	print(('\033[1;32;40m #### Processing scan %d done!')%(int(scan_number[ii])))
end

# Experiment parameters 

# NanoMax convention:
# gamma - horizontal detector
# delta - vertical detector
# gonphi - rotation about vertical axis
# gontheta - rotation about horizontal axis
# photon_energy   = motor_positions.energy
# gonphi          = motor_positions.gonphi # [deg] - can be a range of angles
# gontheta        = gontheta; # [deg] - can be a range of angles

# radius          = motor_positions.radius*1e-3; % [m]
# delta           = -motor_positions.delta;%2.1; % 2.1[deg] these angles are corrected with the sign respecting the rotation rules
# gamma           = motor_positions.gamma; %0.48 12.72[deg] +0.467
# detector_pitch  = 55e-6; % [m]

# nanomax.direct_beam = round([1,  50]);