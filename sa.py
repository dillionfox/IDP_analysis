from sa_calcs import utils
from sa_utils.calc_manager import calc_manager
from sa_utils.data_manager import data_manager
from sa_calcs.sa_prot import sa_prot
from sa_calcs.sa_mem  import sa_mem
from sa_calcs.sa_traj import sa_traj

"""

			################################
			### Structure Analysis Class ###
			###        Dillion Fox       ###
			###          3/2019          ###
			###        UPenn/ORNL        ###
			################################

This class contains some functions that may be useful for analyzing protein structures and membranes.
You can pick and choose the tasks you wish to run by specifying the "calcs" variable upon instantiating 
the class.

For usage, see the examples in the "example_wrappers" folder.

"""

# access previously imported libraries
md = utils.md
np = utils.np
os = utils.os

class SA(sa_prot,sa_traj,sa_mem,calc_manager,data_manager):

	def __init__(self, trajname, top='NULL', name_mod='', outdir='', calcs=[]):
		calc_manager.__init__(self,calcs,outdir,name_mod) 	# make sure the calculations requested make sense and can be done
		data_manager.__init__(self,self.__dict__)	 	# load previously computed data and write data at the end
		sa_traj.__init__(self,trajname,top)			# store trajectory info and some methods that require entire trajectory
		sa_prot.__init__(self)					# main protein attributes and calculations
		sa_mem.__init__(self,outdir,name_mod)			# main membrane attributes and calculations

	def precalcs(self,traj):
		"""
		Calculations that can be done before looping through trajectory

		"""
		self.mem_precalcs(traj,self.precalcs_list)

	def frame_calcs(self,traj,struc,fr):
		"""
		Calculations that are done frame-by-frame

		"""
		if self.protein_analysis:
			prot = struc.atom_slice(struc.topology.select('protein'))[0]
			self.protein_calcs(prot,self.calcs)
		self.diffusion(traj,fr,self.calcs)
		self.membrane_calcs(struc,fr,self.calcs)
		return None

	def traj_calcs(self,traj):
		"""
		Calculations that require all frames at once

		"""
		self.protein_traj_calcs(traj,self.traj_list)
		self.mem_traj_calcs(traj,self.traj_list)
		self.electrostatic_maps(traj,self.traj_list)
		return None

	def post_process(self):
		"""
		All post-processing functions go here

		"""
		self.protein_post_analysis(self.analysis)
		self.traj_post_analysis(self.analysis)
		self.mem_post_analysis(self.analysis)
		return None

	def run(self,mode='default'):
	        """
		Runs and handles all function calls. All data is stored in class objects.

	        """
		from timeit import default_timer as timer
		start = timer()
		global traj_ext ; traj_ext = self.trajname.split('.')[-1]
		if traj_ext not in ['dcd', 'xtc', 'txt']:
			print "do not recognize file type for", self.trajname
			raise ValueError
		#---Print Header
		utils.header()
		#---Code can currently be run in two modes: default, and 'score' mode
		self.mode = mode
		#---Check to see which calculations need to be run
		self.check_input()
		print self.trajname
		#---Load existing data
		if self.mode == 'default' and self.overwrite == 'n': self.load_data()
		elif self.overwrite == 'y': print "OVERWRITING OLD DATA!"
		#---Print log of tasks left to complete
		print "calculations left to do:", self.calcs,self.traj_list
		#---Decide if it's necessary to load trajectories/PDBs
		if len(self.calcs) == 0 and len(self.traj_list) == 0 and len(self.precalcs_list) == 0: LOAD = False
		else: LOAD = True
		#---Run the code
		if self.mode == 'default' and not self.plot_only and LOAD == True:
			print "Loading Trajectory"
			if traj_ext == 'dcd':
				traj = md.load_dcd(self.trajname, self.top)
				self.struc_info(traj[0],len(traj))
			elif traj_ext == 'xtc':
				traj = md.load_xtc(self.trajname, top=self.top)
				self.struc_info(traj[0],len(traj))
			elif traj_ext == 'txt':
				with open(self.trajname) as t:
					self.top = t.readline().rstrip()
					nlines = sum(1 for line in t)
				self.struc_info(md.load(self.top),nlines)
				traj = open(self.trajname)
			#---Only load the necessary frames
			if self.last_frame != -1: traj = traj[self.first_frame:self.last_frame]
			if self.skip_frames != 1: traj = traj[::self.skip_frames]
			#---Pre-calculate some things to avoid computing multiple times
			if traj_ext != 'txt' and len(self.precalcs_list) > 0:
				self.precalcs(traj)
			#---Calculations done on trajectory all at once
			if len(self.traj_list) > 0: self.traj_calcs(traj)
			#---Frame-by-frame calculations
			if len(self.calcs) > 0:
				for fr_,struc in enumerate(traj):
					if self.verbose: print "frame", fr_, "of", len(traj)
					#---If .txt, then the structures have to be loaded one-by-one
					if traj_ext == 'txt': struc = md.load(struc.split("\n")[0])
					#---Run all frame-by-frame calculations
					self.frame_calcs(traj,struc,fr_)
			#---Write data
			self.write_data()
		#---Code can be run in special mode where it references a Rosetta score file and only computes statistics on top N structures
		elif self.mode == 'score': self.run_score_mode() # I copy and pasted the code above without modification. If it doesn't work, move it back.
		#---Run post-processing functions, i.e. plotting, etc.
		self.post_process()
		print "Total execution time:", timer()-start
		return None
