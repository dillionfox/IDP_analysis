from lib_handler import np, md, plt, font, os

"""
FUNCTIONS:

gaussian_chain_2(EED,Rg,outdir,name_mod): None
gaussian_chain(EED): None
semiflexible_chain(EED,outdir,name_mod): None
compute_flory(struc,nres): flory
flory_hist(flory,outdir,name_mod): None

"""

class polymer:
	def __init__(self):
		self.flory = []			# Flory Exponent

	@staticmethod
	def persistence_length(traj):
		pdb = 'temp.pdb'
		traj.save(pdb)
		import MDAnalysis as mda
		from MDAnalysis.analysis import polymer
		
		u = mda.Universe(pdb)
		sel = u.select_atoms('protein and (type C or type N or type O)')

		print "len", len(sel)
		
		p = polymer.PersistenceLength([sel]).run()
		p.perform_fit()
		print "Persistence Length =", p.lp
		os.remove(pdb)

		return None

	@staticmethod
	def persistence_length_2(EED,outdir,name_mod):
		"""
		Compute persistence length
	
		"""
		# Gaussian Chain
		import scipy as sy
		from scipy.optimize import curve_fit
		def EED_dist(r,Nc,p,L):
			#def EED_dist(r,Nc,p):
			return Nc * ((1-np.power(r,2)/(r*L**2)))**(1.5) * ((1-np.power(r,2)/L**2)+1)**(-3) * np.exp(-1.5 * (p*L/(1-np.power(r,2)/L**2)))
		plt.clf()
		n, bins, patches = plt.hist(EED, 50)
		nbins = np.array([(bins[b]+bins[b+1])/2 for b in range(len(bins)-1)])
		p = sy.array([500,8,np.mean(EED)])
		popt, pcov = curve_fit(EED_dist, nbins, n)
		plt.plot(nbins, EED_dist(nbins, *popt), 'r-')
		plt.savefig(outdir+'plength_EED' + name_mod + '.png')
		#np.savetxt(outdir+'plength_EED_fit' + name_mod + '.npy',[nbins,pdf(nbins, *popt)])
		print " !!!!!!!!! Nc,p,L", popt
		return None

	@staticmethod
	def gaussian_chain(EED,Rg,outdir,name_mod):
		"""
		Compute distributions and fit them to Gaussian. If fit is good, it's a Gaussian chain
	
		"""
		from scipy.optimize import curve_fit
		def pdf(rL,R,a):
			return (3.0/(2*np.pi*R**2))**(3.0/2.0) * np.exp(-1.5*((rL-a)/R)**2)
		#---EED
		plt.clf()
		n, bins, patches = plt.hist(EED, 100)
		nbins = [(bins[b]+bins[b+1])/2 for b in range(len(bins)-1)]
		popt, pcov = curve_fit(pdf, nbins, n)
		plt.plot(nbins, pdf(nbins, *popt), 'r-')
		plt.savefig(outdir+'gaussian_chain_EED' + name_mod + '.png')
		np.savetxt(outdir+'gaussian_chain_EED_fit' + name_mod + '.npy',[nbins,pdf(nbins, *popt)])
	
		#--- Rg
		plt.clf()
		n, bins, patches = plt.hist(Rg, 100)
		nbins = [(bins[b]+bins[b+1])/2 for b in range(len(bins)-1)]
		popt, pcov = curve_fit(pdf, nbins, n)
		plt.plot(nbins, pdf(nbins, *popt), 'r-')
		plt.savefig(outdir+'gaussian_chain_Rg' + name_mod + '.png')
		np.savetxt(outdir+'gaussian_chain_Rg_fit' + name_mod + '.npy',[nbins,pdf(nbins, *popt)])
		return None
	
	@staticmethod
	def semiflexible_chain(EED,outdir,name_mod):
		"""
		Compute distributions and compare to skewed Gaussian. If fit is good then it's semi-flexible (like DNA)
	
		"""
	
		# Semiflexible Chain
		from scipy.optimize import curve_fit
		def semiflex(y,Nc,p):
			return Nc * (y)**(-1.5) * (y+1)**(-3) * np.exp(-1.5 * (p/y))
		nbins = 30
		n, bins, patches = plt.hist(EED, 30, normed=1)
		nbins = [(bins[b]+bins[b+1])/2 for b in range(len(bins)-1)]
		nbins = np.asarray(nbins).ravel()
		popt, pcov = curve_fit(semiflex, nbins, n)
		print popt
		plt.plot(nbins, semiflex(nbins, *popt), 'k--')
		plt.xlabel("EED")
		plt.ylabel("P(EED)")
		plt.savefig(outdir+'semiflexible_chain_EED' + name_mod + '.png')
		np.savetxt(outdir+ 'semiflexible_chain_EED_fit' + name_mod + '.npy',[nbins,semiflex(nbins, *popt)])
		return None
	
	def compute_flory(self,struc,nres):
		"""
		The coordinates need to be centered EACH TIME Rg is computed. Therefore you can't just keep adding to S and 
		recomputing the eigenvalues without centering everything. 
	
		"""
		N = range(5,nres-5)
		rg = np.zeros(len(N))
		count = np.zeros(len(N))
		for n in N:
			for r in range(nres-n):
				sel = struc.atom_slice(struc.topology.select('resid ' + str(r) + ' to ' + str(r+n-1)))
				rg[n-5] += md.compute_rg(sel)[0] 
				count[n-5] += 1
		self.flory.append([rg[i]/count[i] for i in range(len(rg))])
		return None
	
	def flory_hist(self,nres,outdir,name_mod):
		"""
		Plot simple histogram of Flory exponents
	
		"""
		from scipy.optimize import curve_fit
		def f(N,b,v):
			return N*v + b
	
		print self.flory
		nframes, npep = np.array(self.flory).shape
		N = np.log(np.asarray(range(5,nres-5)).ravel())
		rg = np.log(np.sum(self.flory,axis=0)/nframes)
	
		plt.clf()
		plt.plot(N,rg,color='b')
		popt, pcov = curve_fit(f, N, rg)
		print popt
	
		err = np.zeros(len(N)) 
		for i in range(len(N)):
			err[i] = np.std(rg[i])
		plt.errorbar(N, rg,yerr=err/(np.sqrt(nframes)-1), fmt='o',color='b')
	
		plt.plot(N, f(N, *popt), 'r-')
		plt.xlabel('log(N)')
		plt.ylabel('log(Rg)')
		plt.text(1.6,-0.14, r'$R_g$ = %.2f*$N^{%.2f}$' % (np.exp(popt[0]),popt[1]), fontdict=font)
		plt.subplots_adjust(top=0.9,bottom=0.15,left=0.18,right=0.85,hspace=0.2,wspace=0.2)
		plt.savefig(outdir+ 'flory_exponents' + name_mod + '.png')
		return None
