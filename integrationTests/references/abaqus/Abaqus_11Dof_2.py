#from job import *
import re
# import numpy as np
from odbAccess import openOdb

# # # START PARAM INPUT

#name of input file
name_inputfile_old='11Dof_1.inp'
name_inputfile_new='11Dof_2.inp' 
name_outputfile='Abaqus_11Dof_2.dat'
name_job='Job2b'

print('script started')

#make lookup table for k
k_body_fl=28e3*0.69
k_tyre_fl=260e3
k_body_fr=28e3*0.69
k_tyre_fr=260e3
k_body_rl=16e3*0.82
k_tyre_rl=260e3
k_body_rr=16e3*0.82
k_tyre_rr=260e3
a=[k_body_fl,k_tyre_fl,k_body_fr,k_tyre_fr,k_body_rl,k_tyre_rl,k_body_rr,k_tyre_rr]
b=0
c=0

l_min = 0.05;
l_max = 0.8;
size_grid = 11;
dl=(l_max-l_min)/(size_grid-1);

L_init=0.3

X=[]
k_grid=[]

for i in range(size_grid):
	X.append(l_min+i*dl)

for i in range(8):
	for j in range(size_grid):
		k_grid.append(a[i]+b*X[j]+c*X[j]*X[j])

with open(name_inputfile_old, 'r') as file:
	nooflines=len(file.readlines())


newfile=''
count=0
idx=0

with open(name_inputfile_old, 'r') as file:	
	for i in range(nooflines):
		line=file.readline()		
		if re.split(',',line)[0] == '*Connector behavior':				
			#find k_index
			if re.split(',',line)[1] == ' name=body_fl-spring-beh\n':
				idx=0
			elif re.split(',',line)[1] == ' name=tyre_fl-spring-beh\n':
				idx=1
			elif re.split(',',line)[1] == ' name=body_fr-spring-beh\n':
				idx=2	
			elif re.split(',',line)[1] == ' name=tyre_fr-spring-beh\n':
				idx=3	
			elif re.split(',',line)[1] == ' name=body_rl-spring-beh\n':
				idx=4	
			elif re.split(',',line)[1] == ' name=tyre_rl-spring-beh\n':
				idx=5	
			elif re.split(',',line)[1] == ' name=body_rr-spring-beh\n':
				idx=6	
			elif re.split(',',line)[1] == ' name=tyre_rr-spring-beh\n':
				idx=7
			else:
				print('unknown behavior name')
			newline = line
			line2 = file.readline() 
			newline = newline + re.split('\n',line2)[0] + ', nonlinear\n'	
			#write k_grid
			for j in range(size_grid):
				newline=newline + '%.6f,%.6f\n' %(k_grid[idx*size_grid+j]*(X[j]-L_init),X[j]-L_init)				
			count=2 #skip the next two lines	
		elif count>1:
			newline=''
			count = count - 1
		else:
			newline=line

		newfile=newfile+newline
	
	
with open(name_inputfile_new, 'w') as file:
	file.write(newfile)	


job=mdb.JobFromInputFile(name_job,name_inputfile_new)
job.submit()
job.waitForCompletion()

odb=openOdb(name_job+'.odb')
inst=odb.rootAssembly.instances['PART-1-1']
step=odb.steps['BE']

noofframes=len(step.frames)
noofnodes=len(inst.nodes)

#history output
rx=step.historyRegions['Node PART-1-1.1'].historyOutputs['UR1'].data
ry=step.historyRegions['Node PART-1-1.1'].historyOutputs['UR2'].data
z=step.historyRegions['Node PART-1-1.1'].historyOutputs['U3'].data

with open(name_outputfile,'w') as file: 		
	file.write('t, z_1, r_x, r_y, z_6, z_7, z_8, z_9, z_10, z_11, z_12, z_13\n')
	frame=step.frames[0]
	for i in range(noofframes):
		frame=step.frames[i]
		t=frame.frameValue
		file.write('%.6f, ' %t)
		file.write('%.12f, %.12f, %.12f, ' %(z[i][1],rx[i][1],ry[i][1]))
		for j in range(noofnodes):
			if j<5: 
				continue
			else:
				val=frame.fieldOutputs['U'].values[j].data[2]
				if j==(noofnodes-1):
					file.write('%.12f ' %val)
				else:
					file.write('%.12f, ' %val)
		file.write('\n')

odb.close()

print('script ended')