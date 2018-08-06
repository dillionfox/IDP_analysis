#!/usr/bin/python
import numpy as np
import sys
import MDAnalysis
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

vec_to_int = {'x':0, 'y':1, 'z':2}
int_to_vec = {0:'x', 1:'y', 2:'z'}
interface_vector = 'z'

def probe(fr,GRO,XTC,cutoff,probe_radius):
	print "frame:", fr, "cutoff:", cutoff
	uni = MDAnalysis.Universe(GRO,XTC)
	uni.trajectory[fr]
	bilayer = uni.select_atoms('resname POPC or resname DOPC or resname DOPS')
	selstring = 'around ' + str(cutoff) +  ' global (name BB)'
	lipids = bilayer.select_atoms(selstring)
	lipid_av = lipids.positions[:,vec_to_int[interface_vector]].mean(axis=0)

	protein_all = uni.select_atoms('protein ')
	resids = np.zeros(len(protein_all.residues.resids))
	selstring = 'around ' + str(cutoff) +  ' global (resname POPC or resname DOPC or resname DOPS)'
	protein = protein_all.select_atoms(selstring)
	protein_av = protein.positions[:,vec_to_int[interface_vector]].mean(axis=0)

	if protein_av > lipid_av: 	# if protein is ABOVE bilayer
		ineq = '<'		# then we are looking for water BELOW protein
		ineq_ = '>'
	else: 				# if protein is BELOW bilayer
		ineq = '>'		# then we are looking for water ABOVE protein
		ineq_ = '<'

	for i in protein.residues.resids:
		pc = protein_all.positions[i-1]
		iv = interface_vector # shorter name
		dl = probe_radius # shorter name
		v0 = vec_to_int[iv]; v1 = (v0+1)%3; v2 = (v0+2)%3
		iv1 = int_to_vec[v1]; iv2 = int_to_vec[v2]

		selstring = '(resname W or name BB) and (prop '+str(iv)+ineq+str(pc[v0])
		selstring += ' and prop ' + str(iv) + ineq_ + str(lipid_av) 
		selstring += ' and prop ' + str(iv1) + ' > ' + str(pc[v1]-dl) 
		selstring += ' and prop ' + str(iv1) + ' < ' + str(pc[v1]+dl) 
		selstring += ' and prop ' + str(iv2) + ' > ' + str(pc[v2]-dl) 
		selstring += ' and prop ' + str(iv2) + ' < ' + str(pc[v2]+dl) 
		selstring += ')'
		water_sel = uni.select_atoms(selstring)

		if len(water_sel) == 0:
			print "CONTACT:", i
			resids[i-1] = 1
	return resids 

def make_map(dat, sequence):
	from matplotlib.ticker import MultipleLocator

	def average_set(curr_list):
		cutoff = 5
	        hist=[]
	        unique_list=sorted(set(curr_list))
	        for el in unique_list:
	                occ=curr_list.count(el)
	                if occ > cutoff:
	                        hist.append(el)
		return hist
	
	fontSize = 10
	res_per_frame = []	
	seq = []
	curr_list = []
	frame_bin = []
	fr_count = 0

	tdat = np.transpose(dat)

	for frame in tdat:
		residues = np.where(frame!=0)[0]
	        for r in residues:
			curr_list.append(r)
	        if fr_count%25==0:
			print "HERE", curr_list
	                hist=average_set(curr_list)
	                for i in range(len(hist)):
	                        frame_bin.append(fr_count)
	                        seq.append(hist[i])
	                res_per_frame.append(len(hist))
	                curr_list=[]
	        fr_count+=1

        fig=plt.figure(1)
        plot=fig.add_subplot(111)
        plt.scatter(frame_bin,seq,s=30,marker='o',lw=0.0)
        plot.tick_params(axis='both', which='minor', labelsize=fontSize)
        #plot.yaxis.set_minor_locator(MultipleLocator(1))
        #plt.grid(b=True,which='major',linestyle='-')
        #plt.grid(b=True,which='minor',linestyle='--',alpha=0.3)
        plot.set_xlabel('Time (ns)',size=fontSize)
        plot.set_ylabel('Residue',size=fontSize)
        plot.set_yticks(np.arange(0,len(sequence)))
        plot.set_yticklabels(sequence,fontSize=fontSize)
        #plt.ylim(26.5,35.5)   
        #plt.xlim(0,900)


def make_hist_by_restype(x, sequence):
	restypes = set(sequence)

	hist = np.zeros(len(restypes))
	resnum = 0
	for i in x:
		resname = sequence[resnum]
		list_index = list(restypes).index(resname)
		hist[list_index]+=sum(i)/sequence.count(resname)
		print resname, sequence.count(resname)
		resnum+=1
	print hist

	hist, restypes = (list(t) for t in zip(*sorted(zip(hist, restypes),reverse=True)))

	fig=plt.figure(1)
	ax=fig.add_subplot(111)
	fig.suptitle('Residue/Lipid Contacts')
	ax.set_xlabel('Residue')
	ax.set_ylabel('Number of Occurrences')

	ax.set_xticks(np.arange(len(restypes)))
	ax.set_xticklabels(restypes,fontSize=5)

	plt.bar(np.arange(len(restypes)),hist,width=1,edgecolor='k')
	plt.xlim(-1,len(restypes))

def make_hist(x, sequence, XTC):
	plt.clf()
	m,n=x.shape
	hist = np.zeros(m)
	resnum = 0
	for i in x:
		hist[resnum]=sum(i)
		resnum+=1

	fig=plt.figure(1)
	ax=fig.add_subplot(111)
	fig.suptitle('Residue/Lipid Contacts')
	ax.set_xlabel('Residue Number')
	ax.set_ylabel('Number of Occurrences')

	ax.set_xticks(np.arange(len(sequence)))
	ax.set_xticklabels(sequence,fontSize=5)

	plt.bar(np.arange(m),hist,width=1,edgecolor='k')
	plt.xlim(-1,len(sequence))
	plt.savefig(XTC.split('.')[0]+"_res_hist.png")

def interface_probe(GRO,XTC,skip_frames,first_frame,last_frame,nthreads,cutoff,probe_radius,sequence):
	print "loading...", XTC
	uni = MDAnalysis.Universe(GRO,XTC)
	print "complete" 
	if last_frame == 'last':
		last_frame = len(uni.trajectory)
	protein = uni.select_atoms('protein')
	nres = len(protein.residues)
	first_resid = protein.residues.resids[0]
	uni = []

	frames = range(first_frame,last_frame,skip_frames)
	resids_per_frame = np.zeros((len(frames),nres))
	contact_frames = np.zeros(len(frames))
	if nthreads == 1:
		resids_per_frame = [probe(fr,GRO,XTC,cutoff,probe_radius) for fr in frames]
	elif nthreads > 1:
		from joblib import Parallel,delayed
		from joblib.pool import has_shareable_memory
		#resids_per_frame = Parallel(n_jobs=nthreads)(delayed(probe,has_shareable_memory)(fr) for fr in frames)
	else:
		print "nthreads can't equal", nthreads

	#print resids_per_frame
	resids_per_frame = np.transpose(resids_per_frame)
	make_hist(resids_per_frame, sequence, XTC)
	#make_map(resids_per_frame, sequence)
	#make_hist_by_restype(resids_per_frame, sequence)

if __name__ == "__main__":

	GRO = sys.argv[1]
	XTC = sys.argv[2]
	skip_frames = 1
	first_frame = 0
	last_frame = 'last'
	nthreads = 12
	cutoff = 5
	probe_radius = 1.0
	sequence = ['R','Q','M','N','F','P','E','R','S','M','D','M','S','N','L','Q','M','D','M','Q','G','R','W','M','D','M','Q','G','R','Y','T','N','P','F','N']

	interface_probe(GRO,XTC,skip_frames,first_frame,last_frame,nthreads,cutoff,probe_radius,sequence)
