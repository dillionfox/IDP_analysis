import numpy as np
import os
import mdtraj as md
import matplotlib as mpl
font = {'family' : 'normal','weight' : 'normal','size'   : 15}
mpl.rc('font', **font)
mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess

#calcs = ["rewrap"]    # rewrap, run, Rg, RMSD, EED
#task = ["PCA"]   # SA calcs
calcs = ["run"]
task = ["PCA","flory","cmaps","SS"]   # SA calcs

xtc_list = "/home/dillion/data/dillion/reflectin/HREX/camp0701/xtc_list.txt" 
tpr_list = "/home/dillion/data/dillion/reflectin/HREX/camp0701/tpr_list.txt" 
camp0701 = [xtc_list, tpr_list]

xtc_list = "/home/dillion/data/dillion/reflectin/HREX/ros00007145-b/xtc_list.txt" 
tpr_list = "/home/dillion/data/dillion/reflectin/HREX/ros00007145-b/tpr_list.txt" 
ros00007145 = [xtc_list, tpr_list]

xtc_list = "/home/dillion/data/dillion/reflectin/HREX/ros00048265-b/xtc_list.txt" 
tpr_list = "/home/dillion/data/dillion/reflectin/HREX/ros00048265-b/tpr_list.txt" 
ros00048265 = [xtc_list, tpr_list]

structures = [camp0701, ros00007145, ros00048265]
#structures = [ros00007145, ros00048265]

global OVERWRITE ; OVERWRITE = False
global skip ; skip = 1000
global start ; start = 0 # ns
global end ; end = 1000     # ns
global first_frame ; first_frame = int(start/2)
global last_frame ; last_frame = int(end/2)
global T; T = [300,313.825,328.286,343.414,359.239,375.794,393.111,411.227,430.177,450]

def rewrap(xtc,tpr,simname):
	if os.path.exists(simname + '_pbc.pdb') and OVERWRITE == False:
		print simname + '_pbc.pdb', "exists"
	elif OVERWRITE == True:
		trjconv = subprocess.Popen(['gmx','trjconv','-s',tpr,'-f',xtc,'-o',simname+'.pdb','-pbc','mol','-ur','rect','-b', str(0), '-e', str(1)],stdin=subprocess.PIPE)
		trjconv.communicate(b'1\n0\n')
		trjconv.wait()

	if os.path.exists(simname + '_pbc.xtc') and OVERWRITE == False:
		print simname + '_pbc.xtc', "exists"
	elif OVERWRITE == True:
		trjconv = subprocess.Popen(['gmx','trjconv','-s',tpr,'-f',xtc,'-o',simname+'.xtc','-pbc','mol','-ur','rect','-skip',str(skip),'-b', str(start*1000), '-e', str(end*1000)],stdin=subprocess.PIPE)
		trjconv.communicate(b'1\n0\n')
		trjconv.wait()

	return None


### Combined PCA ###
def load_PCA_data(name):
	if os.path.isfile(name+'/Rg_3_raw.npy'):
		Rg = np.loadtxt(name+'/Rg_3_raw.npy')
	if os.path.isfile(name+'/EED_3_raw.npy'):
		EED = np.loadtxt(name+'/EED_3_raw.npy')
	if os.path.isfile(name+'/Asph_3_raw.npy'):
		Asph = np.loadtxt(name+'/Asph_3_raw.npy')
	if os.path.isfile(name+'/SASA_3_raw.npy'):
		SASA = np.loadtxt(name+'/SASA_3_raw.npy')
	return [Rg, EED, Asph, SASA]

def run_PCA(Rg,EED,Asph,SASA):
	from IDP_analysis import pca
	outdir = ''
	name_mod = '_combined'
	scores = []
	mode = 'multi'
	trajname = "combined"
	ros_frames = []
	calcs = []
	pca.run_PCA(EED,Rg,SASA,Asph,outdir,name_mod,mode,scores,trajname,ros_frames,calcs)

### RMSD ###
def plot_RMSD(base_name):
	for name in base_name:
		plt.clf()
		for i in range(6):
			j = int(np.floor(i/3))
			run_i = i%3
			RMSD = np.loadtxt(name+"/rmsd_"++str(run_i)+"_pbc.npy")
			plt.plot(RMSD,label=str(run_i))
		plt.legend(loc=2)
		plt.savefig("RMSD_"+name+".png")

def plot_RMSD_panel(base_name):
	plt.clf()
	nrows = 5 ; ncols = 2
	fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, figsize=(6, 6))
	fig.text(0.45, 0.02, 'Time (ns)', ha='center')
	fig.text(0.02, 0.5, 'RMSD (nm)', va='center', rotation='vertical')
	name = base_name[0]
	for i in range(nrows*ncols):
		ax = plt.subplot(nrows,ncols,i+1)
		Rg = np.loadtxt(name+"/rmsd_"+str(i)+"_pbc.npy")
		plt.plot(Rg,label=str(i))

		ax.set_ylim([0,1.9])
		ax.set_xlim([0,251])
		ticks = mpl.ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*2))
		ax.xaxis.set_major_formatter(ticks)

		#if i == ncols - 1:
		#	ax.legend(bbox_to_anchor=(1, 1))
		if i < ncols*(nrows - 1):
			ax.tick_params(axis='x', which='both', bottom='off',top='off',labelbottom='off')
		else:
			ax.tick_params(axis='x', which='both', bottom='on',top='off',labelbottom='on')
		if i%ncols != 0:
			ax.tick_params(axis='y', which='both', left='off',right='off',labelleft='off')
		else:
			ax.tick_params(axis='y', which='both', left='on',right='off',labelleft='on')
		i+=1

	plt.subplots_adjust(top=0.93,bottom=0.1,left=0.1,right=0.78,hspace=0.1,wspace=0.1)
	#plt.savefig("Rg.png")
	plt.show()
	plt.close()
	return None


### Rg ###
def plot_Rg(base_name):
	for name in base_name:
		plt.clf()
		for i in range(6):
			j = int(np.floor(i/3))
			run_i = i%3
			Rg = np.loadtxt(name+"/Rg_"+str(run_i)+"_pbc_raw.npy")
			plt.plot(Rg,label=str(run_i))
		plt.legend(loc=2)
		plt.savefig("Rg_"+name+".png")

def plot_Rg_panel(base_name):
	plt.clf()
	nrows = 5 ; ncols = 2
	fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, figsize=(6, 6))
	fig.text(0.45, 0.02, 'Time (ns)', ha='center')
	fig.text(0.02, 0.5, 'Rg (nm)', va='center', rotation='vertical')
	name = base_name[0]
	for i in range(nrows*ncols):
		ax = plt.subplot(nrows,ncols,i+1)
		Rg = np.loadtxt(name+"/Rg_"+str(i)+"_pbc_raw.npy")[first_frame:last_frame]
		plt.plot(Rg,label=str(i))

		ax.set_ylim([0.8,1.9])
		ax.set_xlim([0,250])
		ticks = mpl.ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*2))
		ax.xaxis.set_major_formatter(ticks)

		#if i == ncols - 1:
		#	ax.legend(bbox_to_anchor=(1, 1))
		if i < ncols*(nrows - 1):
			ax.tick_params(axis='x', which='both', bottom='off',top='off',labelbottom='off')
		else:
			ax.tick_params(axis='x', which='both', bottom='on',top='off',labelbottom='on')
		if i%ncols != 0:
			ax.tick_params(axis='y', which='both', left='off',right='off',labelleft='off')
		else:
			ax.tick_params(axis='y', which='both', left='on',right='off',labelleft='on')
		i+=1

	plt.subplots_adjust(top=0.93,bottom=0.1,left=0.1,right=0.78,hspace=0.1,wspace=0.1)
	plt.savefig(name+"/Rg.png")
	#plt.show()
	plt.close()
	return None

def plot_Rg_hist_panel(base_name):
	plt.clf()
	nrows = 5 ; ncols = 2
	fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, figsize=(6, 6))
	fig.text(0.5, 0.02, 'Rg (nm)', ha='center')
	fig.text(0.02, 0.5, 'Probability', va='center', rotation='vertical')
	name = base_name[0]
	bins = np.linspace(0.8,2.8,40)
	for i in range(nrows*ncols):
		ax = plt.subplot(nrows,ncols,i+1)
		Rg = np.loadtxt(name+"/Rg_"+str(i)+"_pbc_raw.npy")[first_frame:last_frame]
		print len(Rg)
		n, bins, patches = plt.hist(Rg,bins,facecolor='green',alpha=0.75,edgecolor='k',label="T ="+str(int(T[i]))+"K")
		ax.set_ylim([0,60])
		ax.set_xlim([0.9,2.5])

		ax.legend(bbox_to_anchor=(1, 1),frameon=False)
		if i < ncols*(nrows - 1):
			ax.tick_params(axis='x', which='both', bottom='off',top='off',labelbottom='off')
		else:
			ax.tick_params(axis='x', which='both', bottom='on',top='off',labelbottom='on')
		if i%ncols != 0:
			ax.tick_params(axis='y', which='both', left='off',right='off',labelleft='off')
		else:
			ax.tick_params(axis='y', which='both', left='on',right='off',labelleft='on')
		i+=1

	plt.subplots_adjust(top=0.93,bottom=0.1,left=0.1,right=0.9,hspace=0.1,wspace=0.1)
	plt.savefig(name+"/"+name+"_Rg.png")
	#plt.show()
	plt.close()
	return None

### EED ###
def plot_EED_panel(base_name):
	plt.clf()
	nrows = 5 ; ncols = 2
	fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, figsize=(6, 6))
	fig.text(0.45, 0.02, 'Time (ns)', ha='center')
	fig.text(0.02, 0.5, 'EED (nm)', va='center', rotation='vertical')
	name = base_name[0]
	for i in range(nrows*ncols):
		ax = plt.subplot(nrows,ncols,i+1)
		Rg = np.loadtxt(name+"/EED_"+str(i)+"_pbc_raw.npy")[first_frame:last_frame]
		plt.plot(Rg,label=str(i))

		ax.set_ylim([0,9])
		ax.set_xlim([0,250])
		ticks = mpl.ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*2))
		ax.xaxis.set_major_formatter(ticks)

		#if i == ncols - 1:
		#	ax.legend(bbox_to_anchor=(1, 1))
		if i < ncols*(nrows - 1):
			ax.tick_params(axis='x', which='both', bottom='off',top='off',labelbottom='off')
		else:
			ax.tick_params(axis='x', which='both', bottom='on',top='off',labelbottom='on')
		if i%ncols != 0:
			ax.tick_params(axis='y', which='both', left='off',right='off',labelleft='off')
		else:
			ax.tick_params(axis='y', which='both', left='on',right='off',labelleft='on')
		i+=1

	plt.subplots_adjust(top=0.93,bottom=0.1,left=0.1,right=0.78,hspace=0.1,wspace=0.1)
	plt.savefig(name+"/"+name+"_EED.png")
	#plt.show()
	plt.close()
	return None

def plot_EED_hist_panel(base_name):
	plt.clf()
	nrows = 5 ; ncols = 2
	fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, figsize=(6, 6))
	fig.text(0.45, 0.02, 'Time (ns)', ha='center')
	fig.text(0.02, 0.5, 'EED (nm)', va='center', rotation='vertical')
	name = base_name[0]
	bins = np.linspace(0,9,40)
	for i in range(nrows*ncols):
		ax = plt.subplot(nrows,ncols,i+1)
		Rg = np.loadtxt(name+"/EED_"+str(i)+"_pbc_raw.npy")[first_frame:last_frame]
		n, bins, patches = plt.hist(Rg,bins,facecolor='green',alpha=0.75,edgecolor='k',label="T ="+str(int(T[i]))+"K")
		ax.legend(bbox_to_anchor=(1, 1),frameon=False)
		ax.set_ylim([0,35])
		ax.set_xlim([0,10])

		if i < ncols*(nrows - 1):
			ax.tick_params(axis='x', which='both', bottom='off',top='off',labelbottom='off')
		else:
			ax.tick_params(axis='x', which='both', bottom='on',top='off',labelbottom='on')
		if i%ncols != 0:
			ax.tick_params(axis='y', which='both', left='off',right='off',labelleft='off')
		else:
			ax.tick_params(axis='y', which='both', left='on',right='off',labelleft='on')
		i+=1

	plt.subplots_adjust(top=0.93,bottom=0.1,left=0.1,right=0.9,hspace=0.1,wspace=0.1)
	#plt.savefig(name+"/"+name+"_EED.png")
	plt.show()
	plt.close()
	return None

### SS ###
def plot_SS(base_name,H_av,H_err,E_av,E_err,C_av,C_err):
	plt.clf()
	name = base_name[0]
	nrows = 5 ; ncols = 2
	fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, figsize=(6, 6))
	for i in range(nrows*ncols):
		ax = plt.subplot(nrows,ncols,i+1)
		ax.plot(H_av[i],label='Helix',c='r') ;      ax.errorbar(range(len(H_av[i])), H_av[i],yerr=H_err[i], fmt='o',color='r',elinewidth=0.5,markersize=0.01)
		ax.plot(E_av[i],label='Beta Sheet',c='b') ; ax.errorbar(range(len(H_av[i])), E_av[i],yerr=E_err[i], fmt='o',color='b',elinewidth=0.5,markersize=0.01)
		ax.plot(C_av[i],label='Coil',c='k') ;       ax.errorbar(range(len(H_av[i])), C_av[i],yerr=C_err[i], fmt='o',color='k',elinewidth=0.5,markersize=0.01)
		ax.set_ylim([0,1])
		ax.set_xlim([0,35])

		if i == ncols - 1:
			ax.legend(bbox_to_anchor=(1, 1))
		if i < ncols*(nrows - 1):
			ax.tick_params(axis='x', which='both', bottom='off',top='off',labelbottom='off')
		else:
			ax.tick_params(axis='x', which='both', bottom='on',top='off',labelbottom='on')
		if i%ncols != 0:
			ax.tick_params(axis='y', which='both', left='off',right='off',labelleft='off')
		else:
			ax.tick_params(axis='y', which='both', left='on',right='off',labelleft='on')

		i+=1

	#for i in range(nrows*ncols):
	#	plt.text(-130 + (i%4)*48, 8.7 - (np.floor(i/4)%5) * 1.85 ,'(S'+str(i)+')',fontsize=FONT-3)

	plt.subplots_adjust(top=0.93,bottom=0.1,left=0.08,right=0.78,hspace=0.2,wspace=0.2)
	#plt.text(-73, -1.2,'Residue Number',fontsize=FONT)
	#plt.text(-160, 4.85,'Frequency',rotation=90,fontsize=FONT)
	plt.savefig(name+'/SS_'+name+'.png')
	#plt.show()
	plt.close()
	return None

def plot_cmaps(cmaps,name):
	plt.clf()
	nrows = 5 ; ncols = 2
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
	i = 0
	for ax in axes.flat:
		im = ax.imshow(cmaps[i], vmin=0, vmax=1)

		if i < ncols*(nrows - 1):
			ax.tick_params(axis='x', which='both', bottom='off',top='off',labelbottom='off')
		else:
			ax.tick_params(axis='x', which='both', bottom='on',top='off',labelbottom='on')
		if i%ncols != 0:
			ax.tick_params(axis='y', which='both', left='off',right='off',labelleft='off')
		else:
			ax.tick_params(axis='y', which='both', left='on',right='off',labelleft='on')

		i+=1
	#plt.text(-75, -175,'T = ' + str(T) + ' K',fontsize=FONT)
	#plt.text(-75, 55,'Residue Number',fontsize=FONT)
	#plt.text(-160, -90,'Residue Number',rotation=90,fontsize=FONT)
	fig.subplots_adjust(top=0.9,bottom=0.09,left=0.25,right=0.8,hspace=0.15,wspace=0.0)
	cbar_ax = fig.add_axes([0.8, 0.25, 0.01, 0.5])
	cbar = fig.colorbar(im, cax=cbar_ax)

	plt.savefig(name+'/CMAPS_'+name+'.png')
	plt.close()

	return None

def get_simname(dir,simname):
	#return dir+'/traj_'+simname+'_'+str(start)+'_'+str(end)
	return dir+'/traj'+simname

if 'PCA' in calcs:
	Rg_ = [] ; EED_ = [] ; Asph_ = [] ; SASA_ = []

for struc in structures:
	if 'SS' in calcs:
		H = []; He = []; E = []; Ee = []; C = []; Ce = []
	cmaps = []
	base_name = []
	xtc_list, tpr_list = struc
	replica = 0
	for xtc,tpr in zip(open(xtc_list), open(tpr_list)):
		xtc = xtc.split("\n")[0]
		tpr = tpr.split("\n")[0]
		name = xtc.split("/")
		dir = name[-2]
		if not os.path.exists(dir):
			os.mkdir(dir)
		if dir not in base_name:
			base_name.append(dir)
		base = name[-1].split(".")[0].split("p")[1]
		simname = get_simname(dir,base)
		if 'rewrap' in calcs:
			rewrap(xtc,tpr,simname)
		#simname += '_pbc'
		if 'PCA' in calcs and replica == 2:
			Rg, EED, Asph, SASA = load_PCA_data(dir)
			Rg_.append(Rg)
			EED_.append(EED)
			Asph_.append(Asph)
			SASA_.append(SASA)
		if 'SS' in calcs:
			[H_av,H_err] = np.loadtxt(dir+'/'+"SS_H_av_"+str(replica)+".npy")
			[E_av,E_err] = np.loadtxt(dir+'/'+"SS_E_av_"+str(replica)+".npy")
			[C_av,C_err] = np.loadtxt(dir+'/'+"SS_C_av_"+str(replica)+".npy")
			H.append(H_av) ; He.append(H_err)
			E.append(E_av) ; Ee.append(E_err)
			C.append(C_av) ; Ce.append(C_err)
		if 'run' in calcs:
			import sa
			analysis = sa.SA(simname+'.xtc', simname+'.pdb', '_'+base,dir+'/', task)
			analysis.first_frame = first_frame 
			analysis.last_frame =  last_frame
			if OVERWRITE == True:
				analysis.overwrite()
			analysis.run()
		if 'cmaps' in calcs and os.path.isfile(dir+'/'+"CMAPS_"+ str(replica) + "_av.npy") == True:
			##---Contact Maps
			cmaps.append(np.loadtxt(dir+'/'+"CMAPS_"+ str(replica) + "_av.npy"))
		replica += 1
	if 'Rg' in calcs:
		#plot_Rg(base_name)
		#plot_Rg_panel(base_name)
		plot_Rg_hist_panel(base_name)
	if 'RMSD' in calcs:
		#plot_RMSD(base_name)
		plot_RMSD_panel(base_name)
	if 'EED' in calcs:
		#plot_EED_panel(base_name)
		plot_EED_hist_panel(base_name)
	if 'SS' in calcs:
		plot_SS(base_name,H,He,E,Ee,C,Ce)
	if 'cmaps' in calcs:
		plot_cmaps(cmaps, base_name[0])

if 'PCA' in calcs:
	U,V = np.array(Rg_).shape
	N=U*V
	Rg_   = np.array(Rg_).reshape(  N) 
	EED_  = np.array(EED_).reshape( N)
	Asph_ = np.array(Asph_).reshape(N)
	SASA_ = np.array(SASA_).reshape(N)
	run_PCA(Rg_,EED_,Asph_,SASA_)
