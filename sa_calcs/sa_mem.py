from utils import np, md, plt, subprocess, os
from MASA import MASA
from lipid_analysis import lipid_analysis
from mesh_calcs import mesh_calcs
from contacts import contacts

class sa_mem(MASA, lipid_analysis, mesh_calcs, contacts):
	def __init__(self,outdir,name_mod):
		MASA.__init__(self,outdir,name_mod)
		lipid_analysis.__init__(self,outdir,name_mod)
		mesh_calcs.__init__(self,outdir,name_mod)
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
		if 'thickness' in calcs:
			self.thickness(struc,fr)
		if 'interdigitation' in calcs:
			self.interdigit(struc,fr)
		return None

	def mem_precalcs(self,traj,precalcs_list):
		"""
		Things that only need to be computed once at the beginning of a calculation

		'membrane_contacts':self.contacts_precalcs(struc)
		'membrane_analysis':self.make_ndx(self.tpr)
		'av_heights':self.heights_precalcs(struc)

		"""
		print "Running Precalcs", precalcs_list
		if 'membrane_contacts' in precalcs_list:
			self.contacts_precalcs(traj[0])
		if 'av_heights' in precalcs_list or 'thickness' in precalcs_list or 'interdigitation' in precalcs_list:
			self.mesh_precalcs(traj)
			self.protein_mesh(traj)
		if 'av_heights' in precalcs_list:
			self.lipid_mesh(traj,'heights','resname DMPC and not type H')
		if 'thickness' in precalcs_list:
			self.lipid_mesh(traj,'thickness','resname DMPC and name P')
		if 'interdigitation' in precalcs_list:
			self.lipid_mesh(traj,'interdigitation','name "4C31"')
		if 'av_heights' in precalcs_list and 'thickness' in precalcs_list:
			print 'WARNING: THICKNESS MESH WILL OVERWRITE AV_HEIGHTS MESH. RUN CALCS SEPARATELY'
		if 'area_per_lipid' in precalcs_list:
			self.apl_precalcs(traj)
		return None

	def mem_traj_calcs(self,traj,calcs):
		if 'membrane_analysis' in calcs:
			self.make_ndx(self.tpr)
			self.gmx_density(self.trajname,self.tpr,self.overwrite,self.first_frame,self.last_frame)
			self.gmx_order(self.trajname,self.tpr,self.overwrite,self.first_frame,self.last_frame)
		if 'av_interdigitation' in calcs:
			self.compute_av_interdigitation(traj)
		if 'lipid_diffusion' in calcs:
			self.lipid_diffusion(traj)
		if 'bending_modulus' in calcs:
			self.bending_modulus(traj)
		return None

	def mem_post_analysis(self,calcs):
		if 'membrane_contacts' in calcs:
			self.plot_contact_hist(self.outdir,self.name_mod)
		if 'av_heights' in calcs:
			self.normalize_heights()
		if 'thickness' in calcs:
			self.plot_thickness()
		if 'interdigitation' in calcs:
			self.plot_interdigitation()
		if 'membrane_analysis' in calcs:
			self.plot_order()
		if 'av_interdigitation' in calcs:
			self.av_interdigitation_plot()
		if 'membrane_contacts' in calcs and 'av_interdigitation' in calcs:
			self.lipid_phase_trans_plot()
		if 'lipid_diffusion' in calcs:
			self.plot_lipid_diffusion()

	def lipid_phase_trans_plot(self):
		n_contacts = np.sum(self.contact_frames,axis=1)
		plt.clf()
		fig, ax = plt.subplots()
		#cax = ax.scatter(self.av_interdigitation,n_contacts,c=range(len(n_contacts)))
		cax = ax.scatter(self.av_interdigitation,n_contacts,c=range(500))
		fig.colorbar(cax)
		ax.set_xlim([-0.5,0.5])
		ax.set_ylim([-0.5,30])
		plt.savefig(self.outdir+'phase_transition_plot'+self.name_mod+'.png')

