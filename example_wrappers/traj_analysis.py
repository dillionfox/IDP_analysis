import numpy as np
import os
import mdtraj as md
import pytraj as pt
import matplotlib as mpl
font = {'family' : 'normal','weight' : 'normal','size'   : 15}
mpl.rc('font', **font)
#mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess

xtc_list = "/home/dillion/Dropbox/Reflectin/structure_prediction/MD/xtc_list.txt" 
tpr_list = "/home/dillion/Dropbox/Reflectin/structure_prediction/MD/tpr_list.txt" 
gro_list = "/home/dillion/Dropbox/Reflectin/structure_prediction/MD/gro_list.txt" 
calcs = ["run"]
task = ["chain"]
global OVERWRITE ; OVERWRITE = True
global skip ; skip = 1000

ff = {"amber03ws":"_a", "charmm36m":"_c"}

def rewrap(xtc,tpr,simname):
	if os.path.exists(simname + '_pbc.pdb') and OVERWRITE == False:
		print simname + '_pbc.pdb', "exists"
	elif OVERWRITE == True:
		trjconv = subprocess.Popen(['gmx','trjconv','-s',tpr,'-f',xtc,'-o',simname+'_pbc.pdb','-pbc','mol','-ur','rect','-b', str(0), '-e', str(1)],stdin=subprocess.PIPE)
		trjconv.communicate(b'1\n0\n')
		trjconv.wait()

	if os.path.exists(simname + '_pbc.xtc') and OVERWRITE == False:
		print simname + '_pbc.xtc', "exists"
	elif OVERWRITE == True:
		trjconv = subprocess.Popen(['gmx','trjconv','-s',tpr,'-f',xtc,'-o',simname+'_pbc.xtc','-pbc','mol','-ur','rect','-skip',str(skip)],stdin=subprocess.PIPE)
		trjconv.communicate(b'1\n0\n')
		trjconv.wait()

	return None

def shrink(xtc,gro,simname):
	if os.path.exists(simname + '.xtc'):
		print simname + '_ai.xtc', "exists"
	else:
		print "writing", simname + '.xtc'

		traj = md.load_xtc(xtc,top=gro)
		traj = pt.Trajectory(xyz=m_traj.xyz.astype('f8'), top=topology_name)
		traj.autoimage()
		traj.save(simname + '_ai.xtc')

		traj_pt = []

		traj = md.load_xtc(simname + '_ai.xtc',top=gro)
		exit()

		with md.formats.XTCTrajectoryFile(simname + '.xtc', 'w') as f:
			for fr in range(50000)[::100]:
				print "frame", fr
				structure = traj[fr]
				peptide = structure.atom_slice(structure.topology.select('protein'))[0]
				peptide.center_coordinates()
				f.write(peptide.xyz*10)
			peptide.save(simname + '.pdb')
	return None

ff_i = ['c', 'a']
def plot_RMSD(base_name):
	for name in base_name:
		plt.clf()
		for i in range(6):
			j = int(np.floor(i/3))
			run_i = i%3
			RMSD = np.loadtxt(name+"/rmsd_"+ff_i[j]+str(run_i)+"_pbc.npy")
			plt.plot(RMSD,label=ff_i[j]+' '+str(run_i))
		plt.legend(loc=2)
		plt.savefig("RMSD_"+name+".png")

def plot_Rg(base_name):
	for name in base_name:
		plt.clf()
		for i in range(6):
			j = int(np.floor(i/3))
			run_i = i%3
			Rg = np.loadtxt(name+"/Rg_"+ff_i[j]+str(run_i)+"_pbc_raw.npy")
			plt.plot(Rg,label=ff_i[j]+' '+str(run_i))
		plt.legend(loc=2)
		plt.savefig("Rg_"+name+".png")

def plot_Rg_panel(base_name):
	plt.clf()
	nrows = 2 ; ncols = 2
	fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, figsize=(6, 6))
	fig.text(0.45, 0.02, 'Time (ns)', ha='center')
	fig.text(0.02, 0.5, 'Rg (nm)', va='center', rotation='vertical')
	for i in range(nrows*ncols):
		name = base_name[i]
		ax = plt.subplot(nrows,ncols,i+1)
		for n in range(6):
			j = int(np.floor(n/3))
			run_n = n%3
			Rg = np.loadtxt(name+"/Rg_"+ff_i[j]+str(run_n)+"_pbc_raw.npy")
			plt.plot(Rg,label=ff_i[j]+' '+str(run_n))

		ax.set_ylim([0.8,1.9])
		ax.set_xlim([0,50])
		ticks = mpl.ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*2))
		ax.xaxis.set_major_formatter(ticks)

		if i == 1:
			ax.legend(bbox_to_anchor=(1, 1))

		if i < 2:
			ax.tick_params(axis='x', which='both', bottom='off',top='off',labelbottom='off')
		else:
			ax.tick_params(axis='x', which='both', bottom='on',top='off',labelbottom='on')
		if i%2 != 0:
			ax.tick_params(axis='y', which='both', left='off',right='off',labelleft='off')
		else:
			ax.tick_params(axis='y', which='both', left='on',right='off',labelleft='on')
		i+=1

	plt.subplots_adjust(top=0.93,bottom=0.1,left=0.1,right=0.78,hspace=0.1,wspace=0.1)
	#plt.savefig("Rg.png")
	plt.show()
	plt.close()
	return None

def plot_RMSD_panel(base_name):
	plt.clf()
	nrows = 2 ; ncols = 2
	fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, figsize=(6, 6))
	fig.text(0.45, 0.02, 'Time (ns)', ha='center')
	fig.text(0.02, 0.5, 'RMSD (nm)', va='center', rotation='vertical')
	for i in range(nrows*ncols):
		name = base_name[i]
		ax = plt.subplot(nrows,ncols,i+1)
		for n in range(6):
			j = int(np.floor(n/3))
			run_n = n%3
			Rg = np.loadtxt(name+"/rmsd_"+ff_i[j]+str(run_n)+"_pbc.npy")
			plt.plot(Rg,label=ff_i[j]+' '+str(run_n))

		ax.set_ylim([0,1.9])
		ax.set_xlim([0,50])
		ticks = mpl.ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*2))
		ax.xaxis.set_major_formatter(ticks)

		if i == 1:
			ax.legend(bbox_to_anchor=(1, 1))

		if i < 2:
			ax.tick_params(axis='x', which='both', bottom='off',top='off',labelbottom='off')
		else:
			ax.tick_params(axis='x', which='both', bottom='on',top='off',labelbottom='on')
		if i%2 != 0:
			ax.tick_params(axis='y', which='both', left='off',right='off',labelleft='off')
		else:
			ax.tick_params(axis='y', which='both', left='on',right='off',labelleft='on')
		i+=1

	plt.subplots_adjust(top=0.93,bottom=0.1,left=0.1,right=0.78,hspace=0.1,wspace=0.1)
	#plt.savefig("Rg.png")
	plt.show()
	plt.close()
	return None


base_name = []
for xtc,gro,tpr in zip(open(xtc_list), open(gro_list), open(tpr_list)):
	xtc = xtc.split("\n")[0]
	tpr = tpr.split("\n")[0]
	gro = gro.split("\n")[0]
	name = xtc.split("/")
	dir = name[-3]
	if dir not in base_name:
		base_name.append(dir)
	simname = ff[name[-2]]+name[-1].split(".")[0].split("p")[1]
	if 'rewrap' in calcs:
			rewrap(xtc,tpr,dir+'/traj'+simname)
	if 'shrink' in calcs:
		if name[-2] == "charmm36m":
			shrink(xtc,gro,dir+'/traj'+simname)
	simname += '_pbc'
	if 'run' in calcs:
		import sa
		analysis = sa.SA(dir+'/traj'+simname+'.xtc', dir+'/traj'+simname+'.pdb', simname,dir+'/', task)
		analysis.overwrite()
		analysis.run()

if 'Rg' in calcs:
	#plot_Rg(base_name)
	plot_Rg_panel(base_name)
if 'RMSD' in calcs:
	#plot_RMSD(base_name)
	plot_RMSD_panel(base_name)

