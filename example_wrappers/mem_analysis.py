import numpy as np
import os, subprocess
import matplotlib as mpl
font = {'family' : 'normal','weight' : 'normal','size'   : 15}
mpl.rc('font', **font)
if 'DISPLAY' not in os.environ:
	mpl.use('Agg')
import matplotlib.pyplot as plt

TEST = False
tpr_list = "mem_tpr.txt" 
xtc_list = "mem_xtc.txt" 
calcs = ["~rewrap", "~run", "contact_hist"] # rewrap, run
task = ["~PCA","membrane_contacts"]
outdir = "/home/dillion/Dropbox/Reflectin/interface/"
global OVERWRITE ; OVERWRITE = True
global skip ; skip = 1000

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
		if os.path.isfile(outdir+name+'/membrane_contacts_pbc_raw.npy'):
			data_n = np.loadtxt(outdir+name+'/membrane_contacts_pbc_raw.npy')-first_res
			if n % N == 0:
				data = data_n
			else:
				data = np.append(data,data_n[:])
			if (n+1) % N == 0 and n > 1:
				plt.hist(data, bins=range(1,int(max(data))), edgecolor='k', align='left')
				plt.xticks(range(len(seqs[int(np.floor(n/N))])), list(seqs[int(np.floor(n/N))]))
				#plt.show()
				plt.savefig(outdir+name+'/membrane_contacts_all.png')
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

def check_dir(dirname):
	if dirname not in base_name:
		base_name.append(dirname)
	if not os.path.exists(dirname):
		os.makedirs(dirname)
	return None

def detect_ind(name):
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
	if 'rewrap' in calcs: rewrap(xtc,tpr,dir+'/traj'+simname)
	simname += '_pbc'
	if 'run' in calcs:
		import sa
		analysis = sa.SA(dir+'/traj'+simname+'.xtc', dir+'/traj'+simname+'.pdb', simname,outdir+dir+'/', task)
		analysis.overwrite()
		analysis.run()
	if TEST: exit()

if 'contact_hist' in calcs:
	plot_contacts(base_name)
