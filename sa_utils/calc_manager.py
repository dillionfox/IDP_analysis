from calc_lists import calc_lists
import numpy as np
import os

class calc_manager(calc_lists):

	def __init__(self,calcs,outdir,name_mod):
		calc_lists.__init__(self)
		self.verbose = False		# extra print statements
		self.mode = 'default'		# more than one way to run this code: default, cluster, score (rosetta)
		self.overwrite = 'n'
		self.plot_only = False
		self.outdir = outdir		# directory to write output to
		self.name_mod = name_mod	# several functions are reused, so this name differentiates them
		self.calcs = np.array(calcs)	# list of calculations to perform at run time
		self.analysis = []		# list of post-processing functions to be run
		self.traj_list = np.array([])	# list of calcs that require whole trajectory
		self.precalcs_list = np.array([])# list of calculations that run before looping through frames

	def check_calcs(self):
		#---If no calculations are specified, run these
		if len(self.calcs) == 0:
			for c in ['Gyr', 'Rg', 'SASA', 'EED', 'Asph', 'PCA']:
				if c not in self.calcs:
					self.calcs = np.append(self.calcs,c)
		for c in self.calcs:
			#---Check for requests that can't be fulfilled
			if c not in self.known_calcs:
				print c, "is not a known calculation..."
				self.calcs = self.calcs[np.where(self.calcs != c)]
			#---Add dependencies
			if c in self.deps:
				for ext_dep in self.deps[c]:
					if ext_dep not in self.calcs:
						self.calcs = np.append(self.calcs,ext_dep)
			#---Add post-processing funcitons
			if c in self.post_analysis:
				self.analysis.append(c)
			# Pull out functions that don't need to be iterated through individually
			if c in self.traj_master_list:
				self.traj_list = np.append(self.traj_list,c)
				self.calcs = self.calcs[np.where(self.calcs != c)]
			if c in self.precalcs_master_list:
				self.precalcs_list = np.append(self.precalcs_list,c)
		analysis_only = ['chain']
		for c in self.calcs:
			if c in analysis_only:
				self.calcs = self.calcs[np.where(self.calcs != c)]
		return None

	def check_input(self):
		"""
		This function determines which calcs need to be run

		"""
		#---Make sure 'outdir' is formatted properly
		try:
			#---Sometimes I accidentally put / on the end twice
			if self.outdir[-2] == '/' and self.outdir[-1] == '/':
				self.outdir = self.outdir[:-1]
			#---Sometimes I over-correct and don't put any
			if self.outdir[-1] != '/':
				self.outdir = self.outdir+'/'
		except:
			if self.outdir == '':
				pass
			else:
				print "something is weird here. outdir =", self.outdir

		#---If outdir doesn't exist, make it
		if not os.path.exists(self.outdir):
			os.makedirs(self.outdir)

		#---Make sure all calculations are being called properly
		self.check_calcs()

		#---It doesn't make sense to try to overwrite raw data and only plot raw data
		if self.plot_only and self.overwrite:
			print "You can't overwrite raw data and only plot raw data. Continuing, but this will not work."
		#---Diffusion code requires some input
		if 'diffusion' in self.calcs:
			print "reminder: If you're computing the diffusion coefficient from a replica exchange simulation,"
			print "then you must use a continuous trajectory"
			self.R_list = np.array([0, 0.25, 0.5, 0.75, 1.0, 1.25,1.5,1.75])	
			# code is not currently set up to compute diffusion from pdbs. Not hard, just not doing it right now.
			if traj_ext == 'txt':
				print "can't compute diffusion constant from pdb's yet. Change the loop structure to fix it"
				exit()
		#---If ros_frames isn't specified, use 100 by default
		if self.mode == 'score' and self.ros_frames == -1:
			self.ros_frames = 100
		return None
