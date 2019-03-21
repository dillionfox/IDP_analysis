from utils import np, md, plt, subprocess, os
from MASA import MASA
from lipid_analysis import lipid_analysis
from curvature import curvature
from contacts import contacts

class sa_mem(MASA, lipid_analysis, curvature, contacts):
	def __init__(self,outdir,name_mod):
		MASA.__init__(self,outdir,name_mod)
		lipid_analysis.__init__(self,outdir,name_mod)
		curvature.__init__(self,outdir,name_mod)
		contacts.__init__(self,outdir,name_mod)

	def membrane_calcs(self,struc,fr,calcs):
		"""
		Separate function for membrane calcs, which require all atoms, not just protein

		'MASA':self.MASA.append(mem.MASA(struc))
		'membrane_contacts':self.compute_contacts(struc)
		'area_per_lipid':self.compute_area_per_lipid(struc)
		'av_heights':self.heights(struc,fr)

		"""
		if 'MASA' in calcs:
			self.compute_MASA(struc,fr)
		if 'membrane_contacts' in calcs:
			if not self.plot_only:
				self.compute_contacts(struc)
		if 'area_per_lipid' in calcs:
			self.compute_area_per_lipid(struc,fr)
		if 'av_heights' in calcs:
			self.heights(struc,fr)
		return None

	def precalcs(self,traj,precalcs_list):
		"""
		Things that only need to be computed once at the beginning of a calculation

		'membrane_contacts':self.contacts_precalcs(struc)
		'membrane_analysis':self.make_ndx(self.tpr)
		'av_heights':self.heights_precalcs(struc)

		"""
		print "Running Precalcs", precalcs_list
		if 'membrane_contacts' in precalcs_list:
			self.contacts_precalcs(traj[0])
		if 'membrane_analysis' in precalcs_list:
			self.make_ndx(self.tpr)
		if 'av_heights' in precalcs_list:
			self.heights_precalcs(traj[0])
			self.lipid_mesh(traj)
		if 'area_per_lipid' in precalcs_list:
			self.apl_precalcs(traj)
		if 'lipid_phase_transition' in precalcs_list:
			self.lipid_phase_transition(traj[0])
		return None
