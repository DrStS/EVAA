import re
from odbAccess import openOdb

#names
name_inputfile='11Dof_1.inp' 
name_outputfile='Abaqus_11Dof_1.dat'
name_job='Job1'

print('script started')

job=mdb.JobFromInputFile(name_job,name_inputfile)
job.submit()
job.waitForCompletion()

odb=openOdb(name_job + '.odb')
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