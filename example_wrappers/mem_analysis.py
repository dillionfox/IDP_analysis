import numpy as np
import os, subprocess
import matplotlib as mpl
font = {'family' : 'normal','weight' : 'normal','size'   : 15}
mpl.rc('font', **font)
if 'DISPLAY' not in os.environ:
	mpl.use('Agg')
import matplotlib.pyplot as plt

TEST = False
PLOT_ONLY = False
VERBOSE = True
OVERWRITE = False
tpr_list = "mem_tpr.txt" 
xtc_list = "mem_xtc.txt" 
#tpr_list = "new_tpr.txt" 
#xtc_list = "new_xtc.txt" 
calcs = ["~rewrap", "run", "~interdigit","~contact_hist", "~Rg", "~SS"] # rewrap, run
task = ["~MASA","~av_heights", "~membrane_analysis","~membrane_contacts",\
		"~chain","~flory","~area_per_lipid","~PCA","~SS","~Rg",\
		"~interdigitation","emaps"]
outdir = "/home/dillion/Dropbox/Reflectin/interface/"
global skip ; skip = 1000

def plot_SS(base_name):
	"""
	This is super specific for this data set. Will need to modify to accommodate different
	naming schemes.

	"""
	N = 6
	seqs = ['QCNPFNQWSYNRHGCYPGYSYGRNMCY',\
		'QCNPFSQYMNYYGRYWNYPGYNNYYNRNMSYP',\
		'YNSPYWSNWYGRNMYNPCQNNQWYGR']

	data = np.array([])
	for n,name in enumerate(base_name):
		plt.clf()
		if os.path.isfile(outdir+name+'/SS_E_av_pbc_3.5.npy') and os.path.isfile(outdir+name+'/SS_C_av_pbc_3.5.npy') and os.path.isfile(outdir+name+'/SS_H_av_pbc_3.5.npy'):
			E_av, E_err = np.loadtxt(outdir+name+'/SS_E_av_pbc_3.5.npy')
			C_av, C_err = np.loadtxt(outdir+name+'/SS_C_av_pbc_3.5.npy')
			H_av, H_err = np.loadtxt(outdir+name+'/SS_H_av_pbc_3.5.npy')

			if n % N == 0:
				E = E_av
				C = C_av
				H = H_av
			else:
				E += E_av
				C += C_av
				H += H_av
			if (n+1) % N == 0 and n > 1:
				nres = E_av.shape[0]
				plt.plot(H/N,label='Helix',c='r') ;     # plt.errorbar(range(nres), H_av,yerr=H_err/(np.sqrt(500)-1), fmt='o',color='r')
				plt.plot(E/N,label='Beta Sheet',c='b') ;# plt.errorbar(range(nres), E_av,yerr=E_err/(np.sqrt(500)-1), fmt='o',color='b')
				plt.plot(C/N,label='Coil',c='k') ;      # plt.errorbar(range(nres), C_av,yerr=C_err/(np.sqrt(500)-1), fmt='o',color='k')
				plt.legend(loc=1)

				plt.xlabel('Residue')
				plt.ylabel('Probability')
				plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
				plt.savefig(outdir+name+'/membrane_SS_all.png')
				print "plotting", outdir+name+'/membrane_SS_all.png'
	return None

def plot_contacts(base_name):
	"""
	This is super specific for this data set. Will need to modify to accommodate different
	naming schemes.

	"""
	N = 6
	seqs = ['QCNPFNQWSYNRHGCYPGYSYGRNMCY',\
		'QCNPFSQYMNYYGRYWNYPGYNNYYNRNMSYP',\
		'YNSPYWSNWYGRNMYNPCQNNQWYGR']
	first_res = 2

	data = np.array([])
	for n,name in enumerate(base_name):
		plt.clf()
		if os.path.isfile(outdir+name+'/membrane_contacts_pbc_raw.npy'):
			data_n = np.loadtxt(outdir+name+'/membrane_contacts_pbc_raw.npy')-first_res
			if n % N == 0:
				data = data_n
			else:
				data = np.append(data,data_n[:])
			if (n+1) % N == 0 and n > 1:
				#plt.hist(data, bins=range(0,int(max(data))+2), edgecolor='k', align='left')
				hist, bin_edges = np.histogram(data,bins=range(len(seqs[int(np.floor(n/N))])+1))
				plt.bar(bin_edges[:-1], hist/3000.0, width = 1,edgecolor='k')
				plt.xlim(min(bin_edges)-1, max(bin_edges))
				plt.ylim(0,0.8)

				plt.xticks(range(len(seqs[int(np.floor(n/N))])), list(seqs[int(np.floor(n/N))]))
				plt.ylabel('% Frames in Contact')
				plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.1)
				plt.savefig(outdir+name+'/membrane_contacts_all.png')
				print "plotting", outdir+name+'/membrane_contacts_all.png'
	return None

def plot_interdigitation(base_name):
	"""
	This is super specific for this data set. Will need to modify to accommodate different
	naming schemes.

	"""
	N = 6
	seqs = ['QCNPFNQWSYNRHGCYPGYSYGRNMCY',\
		'QCNPFSQYMNYYGRYWNYPGYNNYYNRNMSYP',\
		'YNSPYWSNWYGRNMYNPCQNNQWYGR']

	data = np.array([])
	for n,name in enumerate(base_name):
		fil = outdir+name+'/interdigitation_pbc400-500_raw.npy'
		figname = outdir+name+"/"+fil.split('/')[-1].split('.')[0]+'_all.png'
		if os.path.isfile(fil):
			data_n = np.loadtxt(fil)
			if n % N == 0:
				plt.clf()
			plt.plot(data_n,label=name)
			if (n+1) % N == 0 and n > 1:
				plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
				plt.legend()
				plt.savefig(figname)
				print "plotting", figname
	return None

def plot_Rg(base_name):
	"""
	This is super specific for this data set. Will need to modify to accommodate different
	naming schemes.

	"""
	N = 6
	seqs = ['QCNPFNQWSYNRHGCYPGYSYGRNMCY',\
		'QCNPFSQYMNYYGRYWNYPGYNNYYNRNMSYP',\
		'YNSPYWSNWYGRNMYNPCQNNQWYGR']

	data = np.array([])
	for n,name in enumerate(base_name):
		plt.clf()
		if os.path.isfile(outdir+name+'/Rg_pbc_raw.npy'):
			data_n = np.loadtxt(outdir+name+'/Rg_pbc_raw.npy')
			print "loading", outdir+name+'/Rg_pbc_raw.npy'
			if n % N == 0:
				data = data_n
			else:
				data = np.append(data,data_n[:])
			if (n+1) % N == 0 and n > 1:
				plt.hist(data,bins = np.linspace(0.8,2.8,20), edgecolor='k', align='left')
				plt.xlabel('Rg [nm]')
				plt.xlim(0.5,3.0)
				plt.ylabel('# Occurrences')
				plt.ylim(0,1700)
				plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
				plt.savefig(outdir+name+'/membrane_Rg_all.png')
				print "plotting", outdir+name+'/membrane_Rg_all.png'
	return None

def rewrap(xtc,tpr,simname):

	if os.path.exists(simname + '_pbc.pdb') and OVERWRITE == False:
		print "not making pdb"
		print simname + '_pbc.pdb', "exists"
	else:
		trjconv = subprocess.Popen(['gmx','trjconv','-s',tpr,'-f',xtc,'-o',simname+'_pbc.pdb','-pbc','mol','-ur','rect','-b', str(0), '-e', str(1)],stdin=subprocess.PIPE)
		trjconv.communicate(b'0\n0\n')
		trjconv.wait()

	if os.path.exists(simname + '_pbc.xtc') and OVERWRITE == False:
		print simname + '_pbc.xtc', "exists"
	else:
		trjconv = subprocess.Popen(['gmx','trjconv','-s',tpr,'-f',xtc,'-o',simname+'_pbc.xtc','-pbc','mol','-ur','rect','-dt',str(skip)],stdin=subprocess.PIPE)
		trjconv.communicate(b'0\n0\n')
		trjconv.wait()
	return None

def cat(name,simname):
	print name ; exit()
	for filename in os.listdir(directory):
		if filename.endswith(".xtc") and filename != simname:
			print filename 
		else:
			continue

def check_dir(dirname):
	if dirname not in base_name:
		base_name.append(dirname)
	if not os.path.exists(dirname):
		os.makedirs(dirname)
	return None

def detect_ind(name):
	print name
	if "md.part" in name[-1]:
		return name[-1].split("000")[1].split(".")[0]
	elif "traj_comp" in name[-1] and "000" in name[-1]:
		return name[-1].split(".")[0].split("000")[1]
	elif "traj_comp" in name[-1]:
		return name[-1].split(".")[0].split("p")[-1] 
	else:
		print "make new procedure for creating unique simnames. Using generic for now"
		return str(0)

global base_name ; base_name = []
for xtc,tpr in zip(open(xtc_list), open(tpr_list)):
	xtc = xtc.split("\n")[0] ; tpr = tpr.split("\n")[0]

	# create output folders and determine naming scheme
	name = xtc.split("/")
	check_dir(name[0])
	simname = detect_ind(name)

	# run calculations
	if 'rewrap' in calcs: rewrap(xtc,tpr,name[0]+'/traj'+simname)
	simname += '_pbc'
	if 'cat' in calcs:
		cat(name,simname)
	if 'run' in calcs:
		import sa
		analysis = sa.SA(name[0]+'/traj'+simname+'.xtc', name[0]+'/traj'+simname+'.pdb', simname+'400-500',outdir+name[0]+'/', task)
		analysis.tpr = tpr
		analysis.first_frame = 400
		analysis.last_frame = 500
		#analysis.skip_frames = 100
		if OVERWRITE:
			analysis.overwrite = True
		if PLOT_ONLY:
			analysis.plot_only = True
		if VERBOSE:
			analysis.verbose = True
		analysis.run()
	if TEST: exit()

if 'contact_hist' in calcs:
	plot_contacts(base_name)
if 'Rg' in calcs:
	plot_Rg(base_name)
if 'SS' in calcs:
	plot_SS(base_name)
if 'interdigit' in calcs:
	plot_interdigitation(base_name)
