import os
import subprocess
import numpy as np
import shutil
import random
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory

#---This code can be used several different ways. Choose one or more of the following options.
MODE = ['ANALYZE'] # ABINITIO, EXTRACT, ANALYZE, DEBUG

PATH = '/home/dillion/data/reflectin/structure_prediction/rosetta/N-term/'
FASTA   = PATH+'t000_.fasta.txt'
FRAG3   = PATH+'aat000_03_05.200_v1_3.txt'
FRAG9   = PATH+'aat000_09_05.200_v1_3.txt'
PSIPRED = PATH+'t000_.psipred_ss2.txt'

#---ABINITIO options
global total_structures ; total_structures = 50000
global nruns            ; nruns = 10
global nstrucs          ; nstrucs = int(float(total_structures)/nruns)
global work_dir		; work_dir = os.getcwd()
#---EXTRACT options
global pdb_dir		; pdb_dir = work_dir + "/pdbs/"
#---ANALYZE options
global analysis_dir	; analysis_dir = work_dir + "/analysis/"
global score_file	; score_file = work_dir + "/pdbs/" + "score.fsc"

#---this code relies on a specific scheme for naming directories, defined here
def get_dir(d):
	if d < 10:
		return work_dir+"/0"+str(d)+'/'
	else:
		return work_dir+"/"+str(d)+'/'

#---run Rosetta. If Joblib is available then this will be done in parallel. Otherwise this code is useless.
def run_abinitio(d):
	dir = get_dir(d)
	try:
		os.mkdir(dir)
	except:
		print dir, "already exists!"
	os.chdir(dir)
	shutil.copy2(FASTA, dir+'rob_motif.fasta')
	shutil.copy2(FRAG3, dir+'rob_motif.frag-03_05.200_v1_3')
	shutil.copy2(FRAG9, dir+'rob_motif.frag-09_05.200_v1_3')
	shutil.copy2(PSIPRED, dir+'rob_motif.psipred_ss2')
	rand = random.randint(1,10000000)

	FNULL = open(os.devnull, 'w')
	subprocess.call(['AbinitioRelax', \
		'-in::file::fasta rob_motif.fasta', \
		'-in:file:frag3 rob_motif.frag-03_05.200_v1_3', \
		'-in:file:frag9 rob_motif.frag-09_05.200_v1_3', \
		'-abinitio:relax', \
		'-constant_seed', \
		'-jran ', str(rand), \
		'-relax:fast', \
		'-use_filters true', \
		'-psipred_ss2 rob_motif2_8.psipred_ss2', \
		'-nstruct ', str(nstrucs), \
		'-out:file:silent fast_relax_1.out'], stdout=FNULL, stderr=subprocess.STDOUT)

	return None

#---rename pdbs from each chunk of the results so they range from 0 to N (instead of 0 to N/nruns)
def format_pdbname(i):
	if i < 10:
		return 'refS_0000000'+str(i)+'.pdb'
	elif i < 100:
		return 'refS_000000'+str(i)+'.pdb'
	elif i < 1000:
		return 'refS_00000'+str(i)+'.pdb'
	elif i < 10000:
		return 'refS_0000'+str(i)+'.pdb'
	elif i < 100000:
		return 'refS_000'+str(i)+'.pdb'
	elif i < 1000000:
		return 'refS_00'+str(i)+'.pdb'
	elif i < 10000000:
		return 'refS_0'+str(i)+'.pdb'
	else:
		return str(i)+'.pdb'

#---Run Rosetta's "extract" function
def run_extract(n):
	dir = get_dir(n)
	try:
		os.mkdir(pdb_dir)
	except:
		print "pdbs/ already exists"
	try:
		os.mkdir(dir+"pdbs/")
	except:
		print dir+"pdbs/ already exists"
	FNULL = open(os.devnull, 'w')
	subprocess.call(['extract_pdbs', '-in::file::silent', dir+'fast_relax_1.out', '-out::prefix', dir+'pdbs/ref'], stdout=FNULL, stderr=subprocess.STDOUT)
	for i in range(1,5001):
		j = i+n*5000
		pdbname = format_pdbname(i)
		pdbcopy = format_pdbname(j)
		try:
			shutil.copy2(dir+'pdbs/'+pdbname, pdb_dir + pdbcopy)
		except:
			print dir+'pdbs/'+pdbname, "does not exist"
	return None

#---Each run produces a score file. Combine each, but remove the subsequent headers (keep the first one).
def combine_scorefiles():

	if os.path.exists(score_file):
		print score_file, "exits! Delete and rerun"

	try:
		os.mkdir(pdb_dir)
	except:
		pass

	with open(score_file, "w") as outfile:
		for d in range(0,nruns):
			score_d = get_dir(d) + "score.fsc"
			for it, line in enumerate(open(score_d)):
				if d > 0 and it == 0:
					print "skipping first line"
				outfile.write(line)
	return None

#---make a list of all PDB names
def make_PDB_list():
	pdb_list = open(pdb_dir+"PDB_list.txt", 'w')
	for pdb in sorted(os.listdir(pdb_dir)):
		if pdb.split('.')[1] == 'pdb':
			pdb_list.write(pdb_dir+pdb+'\n')
	pdb_list.close()
	return None

#---run the structure_analysis.py's "sa" class with some task
def analyze(task,mode):
	import sa
	PDB_list = pdb_dir+"PDB_list.txt"
	try:
		os.mkdir(analysis_dir)
	except:
		pass
	analysis = sa.SA(PDB_list,'','',analysis_dir,task)
	analysis.ros_frames = 500
	analysis.score_file = score_file
	if mode == 'score':
		analysis.run(mode)
	else:
		analysis.run()
	return None

if __name__ == "__main__":

	if 'ABINITIO' in MODE:
		runs = Parallel(n_jobs=nruns)(delayed(run_abinitio,has_shareable_memory)(d) for d in range(0,nruns))

	elif 'EXTRACT' in MODE:
		runs = Parallel(n_jobs=nruns)(delayed(run_extract,has_shareable_memory)(d) for d in range(0,nruns))
		combine_scorefiles()
		make_PDB_list()

	elif 'ANALYZE' in MODE:
		analyze(['PCA','cmaps','surface_contacts','flory','chain','SS'],'')

	elif 'DEBUG' in MODE:
		combine_scorefiles()
	else:
		print "keyward", MODE, "not recognized"
