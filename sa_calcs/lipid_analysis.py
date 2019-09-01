from utils import np, md, plt, os, subprocess

class lipid_analysis:
	"""
	This class has several things hard-coded, so watch out.

	Calculate the deuterium order parameter and density of lipids using gmx functions.
	You must already have 'gmx order' and 'gmx density' installed for this to work.

	This code doesn't contribute anything new, it's just a simple wrapper
	so that gmx functions can be added into the pipeline.

	"""

	def __init__(self,outdir,name_mod):
		self.outdir = outdir
		self.name_mod = name_mod
		self.first_frame = -1
		self.last_frame = -1
		self.ndx = self.outdir+'/index_custom.ndx'
		self.ndx_chain = self.outdir+'/chain1.ndx'
		self.order_file = self.outdir+'/order'+self.name_mod+'.xvg'
		self.density_file = self.outdir+'/density'+self.name_mod+'.xvg'
		self.av_interdigitation = None
		self.interdigit_fig = self.outdir+'/av_interdigitation'+self.name_mod+'.png'
		self.diff = None
		self.diff_fig = self.outdir+'/lipid_diffusion'+self.name_mod+'.png'

	def bending_modulus(self,traj):
		"""
		bending modulus

		Not finished...!

		"""
		head_txt = 'resname DMPC and type P'
		head_ind = traj.topology.select(head_txt)
		heads = traj.atom_slice(head_ind)
		head_coors = heads.xyz

		tail_txt = 'name "4C21" or name "4C31"'
		tail_ind = traj.topology.select(tail_txt)
		tails = traj.atom_slice(tail_ind)
		tail_coors = tails.xyz
		tail_coors = np.transpose(np.array([(tail_coors[:,i]+tail_coors[:,i+1])/2.0 for i in range(tail_coors.shape[1])[::2]]),(1,0,2))

		N = (tail_coors-head_coors)[:,:int(len(head_ind)/2),:2] # top leaflet
		grid = np.zeros((9,9))
		

		return None

	def tilt_modulus(self,traj):
		"""
		P(t) = C*sin(t)*exp(-k*t^2/(2*kT))

		a+b^*t^2 = -kT ln(P(t)/sin(t))


		Not finished...!

		"""
		pass

	def lipid_diffusion(self,traj):
		"""
		Calculate planar diffusion of lipid headgroups from Einstein relation
		* Divide by 4*time to get D

		<r^2(t)> = <(r(t)-r(0))^2> = 4Dt

		"""
		sel_ind = traj.topology.select('resname DMPC and type P')
		sel = traj.atom_slice(sel_ind)
		coors = sel.xyz[:,:,:2]
		n_frames, n_lipids, _ = coors.shape
		self.diff = np.zeros(n_frames-1)
		for i in range(n_frames-1):
			self.diff[i] = np.mean(np.linalg.norm(coors[i+1]-coors[i]))
		return None

	def plot_lipid_diffusion(self):
		"""
		Plot diffusion over time

		"""
		plt.clf()
		plt.plot(self.diff)
		print "Plotting", self.diff_fig
		plt.savefig(self.diff_fig)
		return None

	def compute_av_interdigitation(self,traj):
		"""
		This should really be called "interdigitation". The gel phase typically has 
		overlapping (interdigitated) tails.

		"""
		sel_ind = traj.topology.select('name "4C31"')
		sel = traj.atom_slice(sel_ind)
		coors = sel.xyz[:,:,2]
		table, bonds = sel.topology.to_dataframe()
		lipid_ids = [i for i in table['resSeq']]
		n_lipids = len(lipid_ids)
		upper_leaflet = np.mean(coors[:,:n_lipids/2],axis=1)
		lower_leaflet = np.mean(coors[:,n_lipids/2:],axis=1)
		#upper_leaflet = np.min( coors[:,:n_lipids/2],axis=1)
		#lower_leaflet = np.max( coors[:,n_lipids/2:],axis=1)
		self.av_interdigitation = upper_leaflet-lower_leaflet
		return None

	def av_interdigitation_plot(self):
		"""
		Plot interdigitation over time

		"""
		plt.clf()
		plt.plot(self.av_interdigitation)
		print "Plotting", self.interdigit_fig
		plt.savefig(self.interdigit_fig)
		return None

	def make_ndx(self, tpr):
		"""
		Wrapper to make ndx files for 'gmx order' and 'gmx density'. Everything is hard-coded.

		"""
		if not os.path.isfile(self.ndx):
			print "Creating", self.ndx
			print tpr, self.ndx
			trjconv = subprocess.Popen(['gmx','make_ndx','-f',tpr,'-o',self.ndx],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
			trjconv.communicate(b'13 & a C21 | a C22 | a C23 | a C24 | a C25 | a C26 | a C27 | a C28 | a C29 | a C210 | a C211 \
				| a C212 | a C213 | a C214\nname 17 chain1\n13 & a C31 | a C32 | a C33 | a C34 | a C35 | a C36 | a C37 \
				| a C38 | a C39 | a C310 | a C311 | a C312 | a C313 | a C314\nname 18 chain2\nq\n')
			trjconv.wait()
		if not os.path.isfile(self.ndx_chain):
			print "Creating", self.ndx_chain
			trjconv = subprocess.Popen(['gmx','make_ndx','-f',tpr,'-o',self.ndx_chain],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
			trjconv.communicate(b'13 & a C21\n 13 & a C22\n 13 & a C23\n 13 & a C24\n 13 & a C25\n 13 & a C26\n 13 & a C27\n\
				 13 & a C28\n 13 & a C29\n 13 & a C210\n 13 & a C211\n 13 & a C212\n 13 & a C213\n 13 & a C214\n\ndel 0-16\nq\n')
			trjconv.wait()
		return None

	def gmx_density(self,xtc, tpr, overwrite,first_frame=-1,last_frame=-1):
		"""
		Density of lipid atoms as a function of z. Hard-coded atom selection.

		"""
		# gmx density -f traj_pbc.xtc -n index.ndx -center
		if not os.path.isfile(self.density_file):
			if first_frame != -1 and last_frame != -1:
				first_frame *= 1000
				last_frame *= 1000
				trjconv = subprocess.Popen(['gmx','density','-f',xtc,'-n',self.ndx, '-s', tpr, '-center', '-o', self.density_file,'-b',str(first_frame),'-e',str(last_frame)],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
			else:
				trjconv = subprocess.Popen(['gmx','density','-f',xtc,'-n',self.ndx, '-s', tpr, '-center', '-o', self.density_file],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
			trjconv.communicate(b'13\n13\n')
			trjconv.wait()
		return None

	def gmx_order(self,xtc, tpr, overwrite,first_frame=-1,last_frame=-1):
		"""
		Deuterium order parameter for lipids. Hard-coded atom selection.

		"""
		# gmx order -s topol.tpr -f traj_pbc.xtc -n ac2.ndx -d z -od deuter_ac2.xvg
		if not os.path.isfile(self.order_file):
			if first_frame != -1 and last_frame != -1:
				first_frame *= 1000
				last_frame *= 1000
				trjconv = subprocess.Popen(['gmx','order','-f',xtc,'-n',self.outdir+'chain1.ndx','-nr',self.outdir+'chain1.ndx','-s', tpr,'-d','z','-od', self.order_file,'-b',str(first_frame),'-e',str(last_frame)],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
				trjconv.wait()
			else:
				trjconv = subprocess.Popen(['gmx','order','-f',xtc,'-n',self.outdir+'chain1.ndx','-nr',self.outdir+'chain1.ndx','-s', tpr,'-d','z','-od', self.order_file],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
				trjconv.communicate(b'17\n')
				trjconv.wait()
		return None

	def plot_order(self):
		"""
		Another bad name. Plot both 'gmx order' and 'gmx density' data

		"""
		for fil in [self.order_file, self.density_file]:
			x = [] ; y = []
			for line in open(fil):
				try:
					line = line.split()
					line[0]
				except:
					continue
				if line[0] == '#':
					continue
				if len(line) == 2 and '@' not in line[0]:
					x.append(float(line[0]))
					y.append(float(line[1]))

			plt.clf()
			plt.scatter(x,y)
			if fil == self.order_file:
				plt.xlabel('Carbon Number') 
				plt.ylabel('Order Parameter')
				print "plotting", self.outdir+'order_parameter'+self.name_mod+'.png'
				plt.savefig(self.outdir+'order_parameter'+self.name_mod+'.pdf')
			else:
				plt.xlabel('Z') 
				plt.ylabel('Density')
				print "plotting", self.outdir+'lipid_density'+self.name_mod+'.png'
				plt.savefig(self.outdir+'lipid_density'+self.name_mod+'.pdf')
		return None

