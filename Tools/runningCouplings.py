####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module for creating and handling running coupling constants."""

##Import required modules:##
import math
import assertions
import constants
import particleData

print "\n////////////////////////////////"
print "Loading runningCouplings Module:"
print "////////////////////////////////\n"

##Functions:##

def check_is_oneLoopAlphaS(toCheck):
	"""A function to check for an instance of the oneLoopAlphaS class."""
	return isinstance(toCheck,oneLoopAlphaS)

def check_is_oneLoopAlphaEM(toCheck):
	"""A function to check for an instance of the oneLoopAlphaEM class."""
	return isinstance(toCheck,oneLoopAlphaEM)

def check_is_twoLoopAlphaS(toCheck):
	"""A function to check for an instance of the twoLoopAlphaS class."""
	return isinstance(toCheck,twoLoopAlphaS)

def check_is_fourLoopAlphaS(toCheck):
	"""A function to check for an instance of the fourLoopAlphaS class."""
	return isinstance(toCheck,fourLoopAlphaS)

def Nf(QSquared):
	"""A function to return the number of quark flavours that can be produced given a centre of mass energy squared."""
	assert assertions.all_are_numbers([QSquared])
	__QSquared = assertions.force_float_number(QSquared)
	__nf = 0
	for __quarkCode in particleData.knownParticles.get_known_quarks():
		__quarkMass = particleData.knownParticles.get_mass_from_code(__quarkCode)
		if (2.0*__quarkMass < math.sqrt(__QSquared)):
			__nf += 1
	return __nf

def beta0(Nf):
	"""A function to return the value of the first beta function coefficient, beta0."""
	##Source: Paper 'The four-loop QCD beta-function and anomalous dimensions'.
	__Ca, __Tr = constants.Ca(), constants.Tr()
	return ((11.0/3.0)*__Ca - (4.0/3.0)*__Tr*Nf)

def beta1(Nf):
	"""A function to return the value of the second beta function coefficient, beta1."""
	##Source: Paper 'The four-loop QCD beta-function and anomalous dimensions'.
	__Ca, __Tr = constants.Ca(), constants.Tr()
	__Cf, __Ca, __Tr = constants.Cf(), constants.Ca(), constants.Tr()
	return ((34.0/3.0)*__Ca*__Ca - (20.0/3.0)*__Ca*__Tr*Nf - 4.0*__Cf*__Tr*Nf)

def beta2(Nf):
	"""A function to return the value of the third beta function coefficient, beta2."""
	##Source: Paper 'The four-loop QCD beta-function and anomalous dimensions'.
	__Ca, __Tr = constants.Ca(), constants.Tr()
	__Cf, __Ca, __Tr = constants.Cf(), constants.Ca(), constants.Tr()
	__beta2P1 = (2857.0/54.0)*__Ca*__Ca*__Ca
	__beta2P2 = 2.0*__Cf*__Cf*__Tr*Nf - (205.0/9.0)*__Cf*__Ca*__Tr*Nf - (1415.0/27.0)*__Ca*__Ca*__Tr*Nf
	__beta2P3 = (44.0/9.0)*__Cf*__Tr*__Tr*Nf*Nf + (158.0/27.0)*__Ca*__Tr*__Tr*Nf*Nf
	return (__beta2P1 + __beta2P2 + __beta2P3)

def beta3(Nf):
	"""A function to return the value of the fourth beta function coefficient, beta3."""
	##Source: Paper 'The four-loop QCD beta-function and anomalous dimensions'.
	__Ca, __Tr = constants.Ca(), constants.Tr()
	__Nc, __Cf, __Ca, __Tr, __rZ3 = constants.Nc(), constants.Cf(), constants.Ca(), constants.Tr(), constants.aperys_constant()
	__beta3P1 = __Ca*__Ca*__Ca*__Ca * ((150653.0/486.0) - (44.0/9.0)*__rZ3)
	__beta3P2 = __Ca*__Ca*__Ca*__Tr*Nf * ((136.0/3.0)*__rZ3 - (39143.0/81.0))
	__beta3P3 = __Ca*__Ca*__Cf*__Tr*Nf * ((7073.0/243.0) - (656.0/9.0)*__rZ3)
	__beta3P4 = __Ca*__Cf*__Cf*__Tr*Nf * ((352.0/9.0)*__rZ3 - (4204.0/27.0))
	__beta3P5 = 46.0*__Cf*__Cf*__Cf*__Tr*Nf
	__beta3P6 = __Ca*__Ca*__Tr*__Tr*Nf*Nf * ((7930.0/81.0) + (224.0/9.0)*__rZ3)
	__beta3P7 = __Cf*__Cf*__Tr*__Tr*Nf*Nf * ((1352.0/27.0) - (704.0/9.0)*__rZ3)
	__beta3P8 = __Ca*__Cf*__Tr*__Tr*Nf*Nf * ((17152.0/243.0) + (448.0/9.0)*__rZ3)
	__beta3P9 = (424.0/243.0)*__Ca*__Tr*__Tr*__Tr*Nf*Nf*Nf
	__beta3P10 = (1232.0/243.0)*__Cf*__Tr*__Tr*__Tr*Nf*Nf*Nf
	__tensorTermAA = __Nc*__Nc*(__Nc*__Nc + 36.0)/24.0
	__beta3P11 = __tensorTermAA*((704.0/3.0)*__rZ3 - (80.0/9.0))
	__tensorTermFA = __Nc*(__Nc*__Nc + 6.0)/48.0
	__beta3P12 = Nf*__tensorTermFA*((512.0/9.0) - (1664.0/3.0)*__rZ3)
	__denominator1 = 96.0*__Nc*__Nc
	__tensorTermFF = (__Nc*__Nc*__Nc*__Nc - 6.0*__Nc*__Nc + 18) / __denominator1
	__beta3P13 = Nf*Nf*__tensorTermFF*((512.0/3.0)*__rZ3 - (704.0/9.0))
	__beta3H1 = __beta3P1 + __beta3P2 + __beta3P3 + __beta3P4 + __beta3P5 + __beta3P6 + __beta3P7
	__beta3H2 = __beta3P8 + __beta3P9 + __beta3P10 + __beta3P11 + __beta3P12 + __beta3P13
	return (__beta3H1 + __beta3H2)

##Classes:##

class oneLoopAlphaS(object):
	"""A class for handling the one-loop running of the strong coupling constant."""

	def calculate(self,QSquared):
		"""A function for calculating the one-loop strong coupling constant at a given COM energy squared."""
		##Source: 'Determination of the QCD coupling alphaS'.
		assert assertions.all_are_numbers([QSquared])
		__QSquared = assertions.force_float_number(QSquared)
		__Nf= Nf(__QSquared)
		__Mz = particleData.knownParticles.get_mass_from_code(particleData.knownParticles.get_code_from_name('Z-boson'))
		##Encorporate in factor of 1/4Pi to make the beta function match:
		__beta0 = beta0(__Nf) / (4.0*math.pi)
		__denominator1 = (__Mz*__Mz)
		__ln = math.log(__QSquared/__denominator1)
		__denominator2 = 1.0 + (self.__oneLoopAlphaSOfMzSquared * __beta0 * __ln)
		return self.__oneLoopAlphaSOfMzSquared / __denominator2

	def __init__(self):
		"""A function to initiate a one-loop strong coupling constant."""
		self.__oneLoopAlphaSOfMzSquared = constants.one_loop_alphaS_of_Mz_squared()
		self.__showerCutOffEnergy = constants.cut_off_energy()
		self.__alphaSMax = self.calculate(self.__showerCutOffEnergy*self.__showerCutOffEnergy)

	def calculate_weighted(self,QSquared):
		"""A function to return the one-loop strong coupling constant weighted on it's maximum value (at parton shower cutoff)."""
		assert assertions.all_are_numbers([QSquared])
		__QSquared = assertions.force_float_number(QSquared)
		return self.calculate(__QSquared) / self.__alphaSMax

	def get_shower_max(self):
		"""A function to return the maximum one-loop value of alphaS during a parton shower (at the cutoff energy)."""
		return self.__alphaSMax

class oneLoopAlphaEM(object):
	"""A class for handling the one-loop running of the fine structure constant."""

	def calculate(self,QSquared):
		"""A function for calculating the one-loop fine structure constant at a given COM energy squared."""
		##Source: "https://www.ippp.dur.ac.uk/~krauss/Lectures/QuarksLeptons/QCD/AsymptoticFreedom_1.html".
		assert assertions.all_are_numbers([QSquared])
		__QSquared = assertions.force_float_number(QSquared)
		__Mz = particleData.knownParticles.get_mass_from_code(particleData.knownParticles.get_code_from_name('Z-boson'))
		__ln = math.log(__QSquared/(__Mz*__Mz))
		__denominator = 1.0 - (self.__oneLoopAlphaEMOfMzSquared*__ln/(3.0*math.pi))
		return self.__oneLoopAlphaEMOfMzSquared / __denominator

	def __init__(self):
		"""A function to initiate a one-loop fine structure constant."""
		self.__oneLoopAlphaEMOfMzSquared = constants.one_loop_alphaEM_of_Mz_squared()
		self.__showerCutOffEnergy = constants.cut_off_energy()
		##Currently only set for showering on resonance with the Z-boson channel:
		self.__alphaEMMax = self.__oneLoopAlphaEMOfMzSquared

	def calculate_weighted(self,QSquared):
		"""A function to return the one-loop fine structure constant weighted on it's maximum value (at Mz^2)."""
		##Currently only set for showering on resonance with the Z-boson channel.
		assert assertions.all_are_numbers([QSquared])
		__QSquared = assertions.force_float_number(QSquared)
		return self.calculate(__QSquared) / self.__alphaEMMax

	def get_shower_max(self):
		"""A function to return the maximum one-loop value of alphaEM during a parton shower (at Mz^2)."""
		##Currently only set for showering on resonance with the Z-boson channel.
		return self.__alphaEMMax

class twoLoopAlphaS(object):
	"""A class for handling the two-loop running of the strong coupling constant."""

	def calculate(self,Q):
		"""A function for calculating the two-loop strong coupling constant at a given COM energy squared."""
		##Source: Paper 'Measurements of the strong coupling constant and the QCD colour factors using four-jet..
		##..observables from hadronic Z decays', p4, equations 2-5.
		##Unlike the others, these equations are given as a function of Q not Q^2 and so use Mz and alphaS(Mz).
		assert assertions.all_are_numbers([Q])
		__Q = assertions.force_float_number(Q)
		__Nf, __Cf, __twoLoopAlphaSOfMz = Nf(__Q*__Q), constants.Cf(), self.__twoLoopAlphaSOfMz
		__Mz = particleData.knownParticles.get_mass_from_code(particleData.knownParticles.get_code_from_name('Z-boson'))
		##Adjust betas to those of this paper:
		__beta0 = beta0(__Nf) / __Cf
		__beta1 = beta1(__Nf) / (2.0*__Cf*__Cf)
		try:
			assert assertions.safe_division(__Q)
			__ln1 = math.log(__Mz/__Q)
			__wOfQ = 1.0 - (__beta0*__twoLoopAlphaSOfMz*__Cf*__ln1 / (2.0* math.pi))
			__denominator4 = __beta0*2.0*math.pi*__wOfQ
			assert assertions.safe_division(__denominator4)
			__term1 = 1.0 - (__beta1*__twoLoopAlphaSOfMz*__Cf*math.log(__wOfQ)/__denominator4)
			__denominator5 = __wOfQ
			assert assertions.safe_division(__denominator5)
			__result = __twoLoopAlphaSOfMz*__term1 / __denominator5
			return __result
		except:
			##Breaks down below the cut-off energy as expected.
			print "Error in two-loop at Q =", __Q
			return -1.0

	def __init__(self):
		"""A function to initiate a two-loop strong coupling constant."""
		self.__twoLoopAlphaSOfMz = constants.two_loop_alphaS_of_Mz()
		self.__showerCutOffEnergy = constants.cut_off_energy()
		self.__alphaSMax = self.calculate(self.__showerCutOffEnergy*self.__showerCutOffEnergy)

	def calculate_weighted(self,QSquared):
		"""A function to return the two-loop strong coupling constant weighted on it's maximum value (at parton shower cutoff)."""
		assert assertions.all_are_numbers([QSquared])
		__QSquared = assertions.force_float_number(QSquared)
		return self.calculate(__QSquared) / self.__alphaSMax

	def get_shower_max(self):
		"""A function to return the maximum two-loop value of alphaS during a parton shower (at the cutoff energy)."""
		return self.__alphaSMax

class fourLoopAlphaS(object):
	"""A class for handling the four-loop running of the strong coupling constant."""
	##In terms of lambda which has been adjusted so that alphaS(Mz^2) gives 0.118.

	def calculate(self,QSquared):
		"""A function for calculating the four-loop strong coupling constant at a given COM energy squared."""
		##Using approximate analytic solution.
		##Source: 'Quantum Chromodynamics', Dissertori Group, Lpthe Group, P Salam, p3, equation 9.5.
		assert assertions.all_are_numbers([QSquared])
		__QSquared = assertions.force_float_number(QSquared)
		__Nf, __pi = Nf(__QSquared), math.pi
		__beta0 = beta0(__Nf) / (4.0*__pi)
		__beta1 = beta1(__Nf) / ((4.0*__pi)*(4.0*__pi))
		__beta2 = beta2(__Nf) / ((4.0*__pi)*(4.0*__pi)*(4.0*__pi))
		__beta3 = beta3(__Nf) / ((4.0*__pi)*(4.0*__pi)*(4.0*__pi)*(4.0*__pi))
		try:
			__denominator6 = (self.__fourLoopLambda*self.__fourLoopLambda)
			assert assertions.safe_division(__denominator6)
			__ln = math.log(__QSquared/__denominator6)
			__denominator7 = __beta0*__beta0*__ln
			assert assertions.safe_division(__denominator7)
			__term1 = __beta1*math.log(__ln)/__denominator7
			__denominator8 = __beta0*__beta0*__beta0*__beta0*__ln*__ln
			assert assertions.safe_division(__denominator8)
			__term2 = ((__beta1*__beta1 * ((math.log(__ln)*math.log(__ln)) - math.log(__ln) - 1.0)) + __beta0*__beta2) / __denominator8
			__term3P1 = (math.log(__ln)*math.log(__ln)*math.log(__ln)) - (5.0/2.0)*math.log(__ln)*math.log(__ln) - 2.0*math.log(__ln) + 0.5
			__denominator9 = (__beta0*__beta0*__beta0*__beta0*__beta0*__beta0*__ln*__ln*__ln)
			assert assertions.safe_division(__denominator9)
			__term3 = (__beta1*__beta1*__beta1 * __term3P1) / __denominator9
			__term4 = ((3.0*__beta0*__beta1*__beta2*math.log(__ln)) - (0.5*__beta0*__beta0*__beta3)) / __denominator9
			__denominator10 = __beta0*__ln
			assert assertions.safe_division(__denominator10)
			__result = (1.0 - __term1 + __term2 - (__term3 + __term4)) / __denominator10
			return __result
		except:
			print "Error in four-loop at QSquared =", __QSquared
			return -1.0

	def __init__(self):
		"""A function to initiate a four-loop strong coupling constant."""
		self.__fourLoopLambda = constants.four_loop_lambda()
		self.__showerCutOffEnergy = constants.cut_off_energy()
		self.__alphaSMax = self.calculate(self.__showerCutOffEnergy*self.__showerCutOffEnergy)

	def calculate_weighted(self,QSquared):
		"""A function to return the four-loop strong coupling constant weighted on it's maximum value (at parton shower cutoff)."""
		assert assertions.all_are_numbers([QSquared])
		__QSquared = assertions.force_float_number(QSquared)
		return self.calculate(__QSquared) / self.__alphaSMax

	def get_shower_max(self):
		"""A function to return the maximum two-loop value of alphaS during a parton shower (at the cutoff energy)."""
		return self.__alphaSMax

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##
	from matplotlib import pyplot
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes
	from mpl_toolkits.axes_grid1.inset_locator import mark_inset
	import precision

	##Define functions used during testing:##
	def square_function(x):
		assert assertions.all_are_numbers([x])
		"""A function to return the square of a number during testing."""
		return x**2.0

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "////////////////////////////////"
	print "Testing runningCouplings module:"
	print "////////////////////////////////"
	assertions.pause(__name__)
	
	##Setup here:##
	print "\nGenerating test values..."
	tOneLoopAlphaS = oneLoopAlphaS()
	tOneLoopAlphaEM = oneLoopAlphaEM()
	tTwoLoopAlphaS = twoLoopAlphaS()
	tFourLoopAlphaS = fourLoopAlphaS()
	tCutOffEnergy = constants.cut_off_energy()	
	tQRange = range(2,3700,1)
	tNfQRange = range(0,3700000,100)
	tNfs , tOneLoopAlphaSs, tTwoLoopAlphaSs = [], [], []
	tOneLoopAlphaEMs = []
	tFourLoopAlphaSs, tDifferenceInAlphaSs = [], []
	tOneLoopWeightedAlphaSs, tTwoLoopWeightedAlphaSs, tFourLoopWeightedAlphaSs = [], [], []
	tOneLoopWeightedAlphaEMs = []
	tBeta0s, tLns, tDenoms = [] , [], []
	tZBosonMass = particleData.knownParticles.get_mass_from_code(particleData.knownParticles.get_code_from_name('Z-boson'))
	tZBosonMassSquared = tZBosonMass * tZBosonMass
	for tN1, tX1 in enumerate(tQRange):
		tQRange[tN1] = (tQRange[tN1]/10.0)
		tOneLoopAlphaSs.append(tOneLoopAlphaS.calculate(tQRange[tN1]**2))
		tOneLoopAlphaEMs.append(tOneLoopAlphaEM.calculate(tQRange[tN1]**2))
		tTwoLoopAlphaSs.append(tTwoLoopAlphaS.calculate(tQRange[tN1]))
		tFourLoopAlphaSs.append(tFourLoopAlphaS.calculate(tQRange[tN1]**2))
		tDifferenceInAlphaSs.append(tFourLoopAlphaSs[tN1]-tOneLoopAlphaSs[tN1])
		tOneLoopWeightedAlphaSs.append(tOneLoopAlphaS.calculate_weighted(tQRange[tN1]**2))
		tOneLoopWeightedAlphaEMs.append(tOneLoopAlphaEM.calculate_weighted(tQRange[tN1]**2))
		tTwoLoopWeightedAlphaSs.append(tTwoLoopAlphaS.calculate_weighted(tQRange[tN1]))
		tFourLoopWeightedAlphaSs.append(tFourLoopAlphaS.calculate_weighted(tQRange[tN1]**2))
		tBeta0s.append( (11.0 - (2.0/3.0) * Nf(tQRange[tN1]**2))/(4.0*math.pi) )
		tLns.append( math.log( (tQRange[tN1]**2) / (tZBosonMass*tZBosonMass)) )
		tDenoms.append(1.0 + (constants.one_loop_alphaS_of_Mz_squared()*tBeta0s[tN1]*tLns[tN1]))
	print "\nGenerating Nf graph values..."
	tNfQRangeSquared = []
	for tn in range(len(tNfQRange)):
		tNfQRange[tn] = (tNfQRange[tn]/10000.0)
		tNfQRangeSquared.append(tNfQRange[tn]*tNfQRange[tn])
		tNfs.append(Nf(tNfQRange[tn]**2))
	
	tQSquaredRange = map(square_function,tQRange)
	tCutOffEnergySquared = constants.cut_off_energy()*constants.cut_off_energy()

	##Test functions:##
	print "\n--------------------------------------------------\n"
	print "//////////////////"
	print "Testing functions:"
	print "//////////////////"
	assertions.pause(__name__)

	##Test check_is_''LoopAlphaS:##
	print "\n--------------------------------------------------\n"
	print "Testing check_is_''LoopAlphaS:\n"
	tSuccessful = True
	print "Calling check_is_oneLoopAlphaS on instance:" , check_is_oneLoopAlphaS(tOneLoopAlphaS)
	if not check_is_oneLoopAlphaS(tOneLoopAlphaS):
		tSuccessful = False
	print "Calling check_is_oneLoopAlphaS on wrong instance:" , check_is_oneLoopAlphaS(tTwoLoopAlphaS)
	if check_is_oneLoopAlphaS(tTwoLoopAlphaS):
		tSuccessful = False
	print "Calling check_is_oneLoopAlphaS on 1.055:" , check_is_oneLoopAlphaS(1.055)
	if check_is_oneLoopAlphaS(1.055):
		tSuccessful = False
	print "Calling check_is_oneLoopAlphaS on 'word':" , check_is_oneLoopAlphaS('word')
	if check_is_oneLoopAlphaS('word'):
		tSuccessful = False
	print "Calling check_is_oneLoopAlphaEM on instance:" , check_is_oneLoopAlphaEM(tOneLoopAlphaEM)
	if not check_is_oneLoopAlphaEM(tOneLoopAlphaEM):
		tSuccessful = False
	print "Calling check_is_oneLoopAlphaEM on wrong instance:" , check_is_oneLoopAlphaEM(tTwoLoopAlphaS)
	if check_is_oneLoopAlphaEM(tTwoLoopAlphaS):
		tSuccessful = False
	print "Calling check_is_oneLoopAlphaEM on 1.055:" , check_is_oneLoopAlphaEM(1.055)
	if check_is_oneLoopAlphaEM(1.055):
		tSuccessful = False
	print "Calling check_is_oneLoopAlphaEM on 'word':" , check_is_oneLoopAlphaEM('word')
	if check_is_oneLoopAlphaEM('word'):
		tSuccessful = False
	print "Calling check_is_twoLoopAlphaS on instance:" , check_is_twoLoopAlphaS(tTwoLoopAlphaS)
	if not check_is_twoLoopAlphaS(tTwoLoopAlphaS):
		tSuccessful = False
	print "Calling check_is_twoLoopAlphaS on wrong instance:" , check_is_twoLoopAlphaS(tOneLoopAlphaS)
	if check_is_twoLoopAlphaS(tOneLoopAlphaS):
		tSuccessful = False
	print "Calling check_is_twoLoopAlphaS on 1.055:" , check_is_twoLoopAlphaS(1.055)
	if check_is_twoLoopAlphaS(1.055):
		tSuccessful = False
	print "Calling check_is_twoLoopAlphaS on 'word':" , check_is_twoLoopAlphaS('word')
	if check_is_twoLoopAlphaS('word'):
		tSuccessful = False
	print "Calling check_is_fourLoopAlphaS on instance:" , check_is_fourLoopAlphaS(tFourLoopAlphaS)
	if not check_is_fourLoopAlphaS(tFourLoopAlphaS):
		tSuccessful = False
	print "Calling check_is_fourLoopAlphaS on wrong instance:" , check_is_fourLoopAlphaS(tTwoLoopAlphaS)
	if check_is_fourLoopAlphaS(tTwoLoopAlphaS):
		tSuccessful = False
	print "Calling check_is_fourLoopAlphaS on 1.055:" , check_is_fourLoopAlphaS(1.055)
	if check_is_fourLoopAlphaS(1.055):
		tSuccessful = False
	print "Calling check_is_fourLoopAlphaS on 'word':" , check_is_fourLoopAlphaS('word')
	if check_is_fourLoopAlphaS('word'):
		tSuccessful = False
	if tSuccessful:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_is_''LoopAlphaS."
	assertions.pause(__name__)

	##Test Nf function:##
	print "\n--------------------------------------------------\n"
	print "Testing Nf function:\n"
	pyplot.figure()
	pyplot.title(r"$\rm{Number\ of\ available\ quark\ flavours}$",fontsize = 18)
	pyplot.xlabel(r"$\rm{\log{(Q^{2})}}$",fontsize = 16)
	pyplot.ylabel(r"$\rm{N_{f}(Q^{2})}$",fontsize = 16)
	pyplot.plot(tNfQRangeSquared,tNfs,linewidth = 2)
	pyplot.axvline(tCutOffEnergySquared, linestyle = '--', color='red')
	for tn in range(3,7):
		tChangeMass = 2.0 * particleData.knownParticles.get_mass_from_code(tn)
		pyplot.axvline(tChangeMass*tChangeMass, linestyle = '--', color='blue')
		pyplot.text(tChangeMass*tChangeMass*1.1,3.5,r"$\rm{2M_{" + str(particleData.knownParticles.get_name_from_code(tn)[0]) + r"}}$",rotation=90,size = 15)
	pyplot.text(tCutOffEnergySquared,2.7,r"$\rm{Shower\ cut-off}$",rotation=90,size = 15)
	pyplot.xscale('log')
	pyplot.ylim(0.0,6.0)
	pyplot.xlim(0.001,10**6)
	assertions.show_graph()
	print "\nFinished testing Nf function."
	assertions.pause(__name__)

	##Test beta0 function:##
	print "\n--------------------------------------------------\n"
	print "Testing beta0 function:"
	##Using N=3 version from four-loop paper which is the same as Krauss'.
	for tNf1 in range(0,7):
		if (tNf1 == 4):
			assertions.pause(__name__)
		tExpectedBeta0 = 11.0 - (2.0/3.0)*tNf1
		print "\nbeta0("+str(tNf1)+") returns:" , beta0(tNf1)
		print "The expected value is:" , tExpectedBeta0
		if precision.check_numbers_equal(beta0(tNf1),tExpectedBeta0):
			print "These are equal within the code precision: Test successful!"
		else:
			print "Not equal: Test failed!"
	print "\nFinished testing beta0 function."
	assertions.pause(__name__)

	##Test beta1 function:##
	print "\n--------------------------------------------------\n"
	print "Testing beta1 function:"
	##Using N=3 version from four-loop paper which is the same as Krauss'.
	for tNf2 in range(0,7):
		if (tNf2 == 4):
			assertions.pause(__name__)
		tExpectedBeta1 = 102.0 - (38.0/3.0)*tNf2
		print "\nbeta1("+str(tNf2)+") returns:" , beta1(tNf2)
		print "The expected value is:" , tExpectedBeta1
		if precision.check_numbers_equal(beta1(tNf2),tExpectedBeta1):
			print "These are equal within the code precision: Test successful!"
		else:
			print "Not equal: Test failed!"
	print "\nFinished testing beta1 function."
	assertions.pause(__name__)

	##Test beta2 function:##
	print "\n--------------------------------------------------\n"
	print "Testing beta2 function:"
	##Using N=3 version from four-loop paper which is the same as Krauss'.
	for tNf3 in range(0,7):
		if (tNf3 == 4):
			assertions.pause(__name__)
		tExpectedBeta2 = (2857.0/2.0) - (5033.0/18.0)*tNf3 + (325.0/54.0)*tNf3*tNf3
		print "\nbeta2("+str(tNf3)+") returns:" , beta2(tNf3)
		print "The expected value is:" , tExpectedBeta2
		if precision.check_numbers_equal(beta2(tNf3),tExpectedBeta2):
			print "These are equal within the code precision: Test successful!"
		else:
			print "Not equal: Test failed!"
	print "\nFinished testing beta2 function."
	assertions.pause(__name__)

	##Test beta3 function:##
	print "\n--------------------------------------------------\n"
	print "Testing beta3 function:"
	##Using N=3 version from four-loop paper which is presumed the same as Krauss'.
	for tNf4 in range(0,7):
		if (tNf4 == 4):
			assertions.pause(__name__)
		tExpectedBeta3P1 = (149753.0/6.0) + 3564.0*constants.aperys_constant()
		tExpectedBeta3P2 = ((1078361.0/162.0) + (6508.0/27.0)*constants.aperys_constant())*tNf4
		tExpectedBeta3P3 = ((50065.0/162.0) + (6472.0/81.0)*constants.aperys_constant())*(tNf4**2)
		tExpectedBeta3P4 = (1093.0/729.0)*(tNf4**3)
		tExpectedBeta3 = tExpectedBeta3P1 - tExpectedBeta3P2 + tExpectedBeta3P3 + tExpectedBeta3P4
		print "\nbeta3("+str(tNf4)+") returns:" , beta3(tNf4)
		##Beta3 should always be positive for positive nf according to the four-loop paper.
		print "For nf > 0, check beta3 positive:" , (beta3(tNf4) >= 0)
		print "The expected value is:" , tExpectedBeta3
		if precision.check_numbers_equal(beta3(tNf4),tExpectedBeta3):
			print "These are equal within the code precision: Test successful!"
		else:
			print "Not equal: Test failed!"
	print "\nFinished testing beta3 function."
	assertions.pause(__name__)

	##Test alphaS class':##
	print "\n--------------------------------------------------\n"
	print "//////////////////////"
	print "Testing alpha classes:"
	print "//////////////////////"
	assertions.pause(__name__)

	##__init__ functions tested implicitly.

	##Test get_shower_max() functions:##
	print "\n--------------------------------------------------\n"
	print "Testing get_shower_max() functions:\n"
	tSuccessful = True
	print "One loop gives:", tOneLoopAlphaS.get_shower_max()
	print "and we expect:", tOneLoopAlphaS.calculate(tCutOffEnergy*tCutOffEnergy)
	if (not precision.check_numbers_equal(tOneLoopAlphaS.get_shower_max(),tOneLoopAlphaS.calculate(tCutOffEnergy*tCutOffEnergy))):
		tSuccessful = False
	print "One loop EM gives:", tOneLoopAlphaEM.get_shower_max()
	print "and we expect:", tOneLoopAlphaEM.calculate(tZBosonMassSquared)
	if (not precision.check_numbers_equal(tOneLoopAlphaEM.get_shower_max(),tOneLoopAlphaEM.calculate(tZBosonMassSquared))):
		tSuccessful = False
	print "Two loop gives:", tTwoLoopAlphaS.get_shower_max()
	print "and we expect:", tTwoLoopAlphaS.calculate(tCutOffEnergy*tCutOffEnergy)
	if (not precision.check_numbers_equal(tTwoLoopAlphaS.get_shower_max(),tTwoLoopAlphaS.calculate(tCutOffEnergy*tCutOffEnergy))):
		tSuccessful = False
	print "Four loop gives:", tFourLoopAlphaS.get_shower_max()
	print "and we expect:", tFourLoopAlphaS.calculate(tCutOffEnergy*tCutOffEnergy)
	if (not precision.check_numbers_equal(tFourLoopAlphaS.get_shower_max(),tFourLoopAlphaS.calculate(tCutOffEnergy*tCutOffEnergy))):
		tSuccessful = False
	if tSuccessful:
		print "\nAll are equal within the code precision: Test successful!"
	else:
		print "\nNOT all equal within the code precision: Test Failed!"
	print "\nFinished testing get_shower_max() functions."
	assertions.pause(__name__)

	##Test one-loop jump direction:##
	print "\n--------------------------------------------------\n"
	print "Testing one-loop jump direction:\n"
	pyplot.figure()
	pyplot.plot(tQRange,tBeta0s, label = r"$\rm{\beta_{0}}$",linewidth = 2)
	pyplot.plot(tQRange,tLns, label = r"$\rm{\ln(Q^{2}\ /\ M_{z}^{2})}$",linewidth = 2)
	pyplot.plot(tQRange,tDenoms, label = r"$\rm{Denominator}$",linewidth = 2)
	pyplot.axvline(tCutOffEnergySquared, linestyle = '--', color = 'red')
	pyplot.text(tCutOffEnergySquared,-7, r"$\rm{Shower\ cut-off}$",rotation=90, size = 14)
	pyplot.title(r"$\rm{One-loop\ terms\ for\ the\ \alpha_{s}}$",fontsize = 18)
	pyplot.xlabel(r"$\rm{(Q)}$",fontsize = 16)
	pyplot.ylabel(r"$\rm{F(Q)}$",fontsize = 16)
	pyplot.legend(loc = "lower right")
	pyplot.xlim(-10,370)
	assertions.show_graph()
	print "\nFinished testing one-loop jump direction."
	assertions.pause(__name__)

	##Test one-loop calculate function:##
	print "\n--------------------------------------------------\n"
	print "Testing one-loop calculate function:\n"
	tFig, tAx = pyplot.subplots() ##Create a new figure with a default 111 subplot.
	tAx.plot(tQSquaredRange,tOneLoopAlphaSs,linewidth = 2)
	tAxins = inset_axes(tAx, 2, 2, loc=1) ##Zoom-factor: 2.5, location: upper-left.
	tAxins.plot(tQSquaredRange,tOneLoopAlphaSs,linewidth = 2)
	tX1, tX2, tY1, tY2 = 4, 12, 0.2, 0.4 ##Specify the limits.
	tAxins.set_xlim(tX1, tX2) ##Apply the x-limits.
	tAxins.set_ylim(tY1, tY2) ##Apply the y-limits.
	tAxins.set_xscale('log')
	tAx.set_xscale('log')
	tAx.set_ylim(0.0,3.6)
	tAx.set_xlim(0.1,15600)
	tAx.set_title(r"$\rm{One-lo}op\ \alpha_{s}$",fontsize = 18)
	tAx.set_xlabel(r"$\log(Q^{2})$",fontsize = 16)
	tAx.set_ylabel(r"$\rm{\alpha_{S}\ (Q^{2})}$",fontsize = 16)
	tAx.axvline(tCutOffEnergySquared, linestyle = '--', color = 'red')
	tAx.text(tCutOffEnergySquared+0.1,3, r"$\rm{Shower\ cut-off}$",rotation=90, size = 14)
	pyplot.yticks(visible=False)
	pyplot.xticks(visible=False)
	mark_inset(tAx, tAxins, loc1=2, loc2=4, fc="none", ec="0.5")
	assertions.show_graph()
	print "\nFinished testing one-loop calculate function."
	assertions.pause(__name__)

	##Test EM one-loop calculate function:##
	print "\n--------------------------------------------------\n"
	print "Testing EM one-loop calculate function:\n"
	tFig, tAx = pyplot.subplots() ##Create a new figure with a default 111 subplot.
	tAx.plot(tQSquaredRange,tOneLoopAlphaEMs,linewidth = 2)
	tAx.set_ylim(0.00768,0.00778)
	tAx.set_xlim(-2000,138000)
	tAx.set_title(r"$\rm{One-loop\ \alpha_{EM}}$",fontsize = 18)
	tAx.set_xlabel(r"$Q^{2}$",fontsize = 16)
	tAx.set_ylabel(r"$\rm{\alpha_{EM}\ (Q^{2})}$",fontsize = 16)
	tAx.axvline(tCutOffEnergySquared, linestyle = '--', color = 'red')
	tAx.text(tCutOffEnergySquared+1000,0.00772, r"$\rm{Shower\ cut-off}$",rotation=90, size = 14)
	assertions.show_graph()
	print "\nFinished testing EM one-loop calculate function."
	assertions.pause(__name__)

	##Test two-loop calculate function:##
	print "\n--------------------------------------------------\n"
	print "Testing two-loop calculate function:\n"
	pyplot.figure()
	pyplot.plot(tQRange,tTwoLoopAlphaSs,linewidth = 2)
	pyplot.axvline(tCutOffEnergySquared, linestyle = '--', color = 'red')
	pyplot.ylim([-0.25,1.0])
	pyplot.text(tCutOffEnergySquared,1.9, r"$\rm{Shower\ cut-off}$",rotation=90, size = 14)
	pyplot.title(r"$\rm{Two-loop\ \alpha_{s}}$",fontsize = 18)
	pyplot.xlabel(r"$\rm{\log(Q)}$",fontsize = 16)
	pyplot.ylabel(r"$\rm{\alpha_{S}\ (Q^{2})}$",fontsize = 16)
	pyplot.xscale('log')
	pyplot.ylim(0.0,2.0)
	assertions.show_graph()
	print "\nFinished testing two-loop calculate function."
	assertions.pause(__name__)

	##Test four-loop calculate function:##
	print "\n--------------------------------------------------\n"
	print "Testing four-loop calculate function:\n"
	pyplot.figure()
	pyplot.plot(tQRange,tFourLoopAlphaSs,linewidth = 2)
	pyplot.axvline(tCutOffEnergySquared, linestyle = '--', color = 'red')
	pyplot.ylim([-0.25,1.0])
	pyplot.text(tCutOffEnergySquared,8, r"$\rm{Shower\ cut-off}$",rotation=90, size = 14)
	pyplot.title(r"$\rm{Four-loop\ \alpha_{s}}$",fontsize = 18)
	pyplot.xlabel(r"$\rm{\log(Q)}$",fontsize = 16)
	pyplot.ylabel(r"$\rm{\alpha_{S}\ (Q^{2})}$",fontsize = 16)
	pyplot.xscale('log')
	pyplot.ylim(0.0,10.5)
	assertions.show_graph()
	print "\nFinished testing four-loop calculate function."
	assertions.pause(__name__)

	##Compare all:##
	print "\n--------------------------------------------------\n"
	print "Comparing all:\n"
	pyplot.figure()
	pyplot.plot(tQRange,tOneLoopAlphaSs,label=r"$\rm{One-loop}$",linewidth = 2)
	pyplot.plot(tQRange,tOneLoopAlphaEMs,label=r"$\rm{One-loop\ EM}$",linewidth = 2)
	pyplot.plot(tQRange,tTwoLoopAlphaSs,label=r"$\rm{Two-loop}$",linewidth = 2)
	pyplot.plot(tQRange,tFourLoopAlphaSs,label=r"$\rm{Four-loop}$",linewidth = 2)
	pyplot.axvline(tCutOffEnergySquared, linestyle = '--', color = 'red')
	pyplot.ylim([-0.25,1.0])
	pyplot.text(tCutOffEnergySquared,3.5, r"$\rm{Shower\ cut-off}$",rotation=90, size = 14)
	pyplot.title(r"$\rm{Comparing\ \alpha_{s}\ at\ different\ loop\ orders}$",fontsize = 18)
	pyplot.xlabel(r"$\rm{log(Q)}$",fontsize = 16)
	pyplot.ylabel(r"$\rm{\alpha_{S}\ (Q^{2})}$",fontsize = 16)
	pyplot.legend(loc=1)
	pyplot.xscale('log')
	pyplot.ylim(0.0,4.0)
	assertions.show_graph()
	print "\nFinished comparing all."
	assertions.pause(__name__)

	##Compare four-loop and one-loop function:##
	print "\n--------------------------------------------------\n"
	print "Comparing four-loop and one-loop function:\n"
	pyplot.figure()
	pyplot.plot(tQRange,tDifferenceInAlphaSs,linewidth = 2)
	pyplot.axvline(tCutOffEnergySquared, linestyle = '--', color = 'red')
	pyplot.ylim([-0.25,1.0])
	pyplot.text(tCutOffEnergySquared,0.3, r"$\rm{Shower\ cut-off}$",rotation=90, size = 14)
	pyplot.title(r"$\rm{Four-loop\ correction\ to\ the\ one-loop\ \alpha_{s}}$",fontsize = 18)
	pyplot.xlabel(r"$\rm{\log(Q)}$",fontsize = 16)
	pyplot.ylabel(r"$\rm{Difference\ in\ \alpha_{S}\ (Q^{2})}$",fontsize = 16)
	pyplot.xscale('log')
	pyplot.ylim(-0.2,0.4)
	assertions.show_graph()
	print "\nFinished comparing four-loop and one-loop function."
	assertions.pause(__name__)

	##Test all weighted_''_loop functions:##
	print "\n--------------------------------------------------\n"
	print "Testing all weighted_''_loop functions:\n"
	pyplot.figure()
	pyplot.plot(tQRange,tOneLoopWeightedAlphaSs,label=r"$\rm{One-loop}$",linewidth = 2)
	pyplot.plot(tQRange,tOneLoopWeightedAlphaEMs,label=r"$\rm{One-loop\ EM}$",linewidth = 2)
	pyplot.plot(tQRange,tTwoLoopWeightedAlphaSs,label=r"$\rm{Two-loop}$",linewidth = 2)
	pyplot.plot(tQRange,tFourLoopWeightedAlphaSs,label=r"$\rm{Four-loop}$",linewidth = 2)
	pyplot.axvline(tCutOffEnergySquared, linestyle = '--', color = 'red')
	pyplot.ylim([-0.25,1.0])
	pyplot.text(tCutOffEnergySquared,1.7,r"$\rm{Shower\ cut-off}$",rotation=90, size = 14)
	pyplot.title(r"$\rm{Weighted\ \alpha_{S/EM}\ at\ different\ loop\ orders}$",fontsize = 18)
	pyplot.xlabel(r"$\rm{\log(Q)}$",fontsize = 16)
	pyplot.ylabel(r"$\rm{\alpha_{S/EM}\ (Q^{2})\ /\ \alpha_{S/EM}\ (Q^{2})_{max}}$",fontsize = 16)
	pyplot.legend(loc=1)
	pyplot.xscale('log')
	pyplot.ylim(0.0,2.0)
	assertions.show_graph()
	print "\nFinished testing all weighted_''_loop functions."
	assertions.pause(__name__)

	##Generating alphaS report graph:##
	print "\n--------------------------------------------------\n"
	print "Generating alphaS report graph:\n"
	pyplot.figure()
	pyplot.plot(tQSquaredRange,tOneLoopWeightedAlphaSs,label=r'$\rm{One-loop}$',linewidth = 2)
	pyplot.axvline(tCutOffEnergySquared, linestyle = '--', color = 'red')
	pyplot.ylim([0.0,1.2])
	pyplot.xlim([0.7,tQSquaredRange[-1]])
	pyplot.text(tCutOffEnergySquared,0.6, r"$\rm{Shower\ cut-off}$",rotation=90, fontsize = 18)
	#pyplot.title(r"$\rm{One\ Loop\ Weighted\ \alpha_{s}}$",fontsize = 24)
	pyplot.xlabel(r"$\rm{\log(Q^{2})}$",fontsize = 22)
	pyplot.ylabel(r"$\rm{\alpha_{s}(Q^{2})/\alpha_{s}(Q_{max}^{2})}$",fontsize = 22)
	pyplot.xscale('log')
	pyplot.xticks(fontsize = 15)
	pyplot.yticks(fontsize = 15)
	assertions.show_graph()
	print "\nFinished Generating alphaS report graph."
	assertions.pause(__name__)

	##Create alphaS seminar graph:##
	print "\n--------------------------------------------------\n"
	print "Creating alphaS seminar graph:\n"
	tFig, tAx = pyplot.subplots() ##Create a new figure with a default 111 subplot.
	tAx.plot(tQSquaredRange,tOneLoopWeightedAlphaSs,linewidth = 2)
	tAxins = inset_axes(tAx, 2, 3, loc=1) ##Zoom-factor: 2.5, location: upper-left.
	tAxins.plot(tQSquaredRange,tOneLoopWeightedAlphaSs,linewidth = 2)
	tX1, tX2, tY1, tY2 = 1, 150, 0.3, 0.7 ##Specify the limits.
	tAxins.set_xlim(tX1, tX2) ##Apply the x-limits.
	tAxins.set_ylim(tY1, tY2) ##Apply the y-limits.
	tAxins.set_xscale('log')
	tAx.set_xscale('log')
	tAx.set_ylim(0.0,1.2)
	tAx.set_xlim(0.2,tQSquaredRange[-1])
	tAx.set_title(r"$\rm{Weighted\ one-loop\ \alpha_{s}}$",fontsize = 18)
	tAx.set_xlabel(r"$\log(Q^{2})$",fontsize = 16)
	tAx.set_ylabel(r"$\rm{\alpha_{s}(Q^{2})/\alpha_{s}(Q_{max}^{2})}$",fontsize = 16)
	tAx.axvline(tCutOffEnergySquared, linestyle = '--', color = 'red')
	tAx.text(0.6,0.4, r"$\rm{Shower\ cut-off}$",rotation=90, size = 14)
	pyplot.yticks(visible=False)
	pyplot.xticks(visible=False)
	mark_inset(tAx, tAxins, loc1=2, loc2=4, fc="none", ec="0.5")
	assertions.show_graph()
	print "\nFinished Creating alphaS seminar graph."
	assertions.pause(__name__)

	##Test reference values:##
	print "\n--------------------------------------------------\n"
	print "Testing reference values:\n"
	print "Z-boson mass is:", tZBosonMass
	print "The one-loop value of alphS(Mz^2) used was:" , constants.one_loop_alphaS_of_Mz_squared()
	print "   Calculating one-loop alphaS(Mz^2) gives:" , tOneLoopAlphaS.calculate(tZBosonMassSquared)
	print "The one-loop value of alphEM(Mz^2) used was:" , constants.one_loop_alphaEM_of_Mz_squared()
	print "   Calculating one-loop alphaEM(Mz^2) gives:" , tOneLoopAlphaEM.calculate(tZBosonMassSquared)
	print "The two-loop value of alphaS(Mz) used was:" , constants.two_loop_alphaS_of_Mz()
	print "   Calculating two-loop alphaS(Mz) gives:" , tTwoLoopAlphaS.calculate(tZBosonMass)
	print "The four-loop value of lambda used was:" , constants.four_loop_lambda()
	print "   Calculating four-loop alphaS(Mz^2) gives:" , tFourLoopAlphaS.calculate(tZBosonMassSquared)
	if precision.check_numbers_equal(constants.one_loop_alphaS_of_Mz_squared(),tOneLoopAlphaS.calculate(tZBosonMassSquared)):
		print "\nOne-loop equal within the code precision: Test successful!"
	else:
		print "\nOne-loop NOT equal within the code precision: Test Failed!"
	if precision.check_numbers_equal(constants.one_loop_alphaEM_of_Mz_squared(),tOneLoopAlphaEM.calculate(tZBosonMassSquared)):
		print "\nOne-loop EM equal within the code precision: Test successful!"
	else:
		print "\nOne-loop EM NOT equal within the code precision: Test Failed!"
	if precision.check_numbers_equal(constants.two_loop_alphaS_of_Mz(),tTwoLoopAlphaS.calculate(tZBosonMass)):
		print "\nTwo-loop equal within the code precision: Test successful!"
	else:
		print "\nTwo-loop NOT equal within the code precision: Test Failed!"
	if precision.check_numbers_equal(constants.one_loop_alphaS_of_Mz_squared(),tFourLoopAlphaS.calculate(tZBosonMassSquared)):
		print "\nFour-loop equal within the code precision: Test successful!"
	else:
		print "\nFour-loop NOT equal within the code precision: Test Failed!"
	print "\nFinished testing reference values."
	assertions.pause(__name__)

	##Done testing:##
	print"\n---------------------------------------------\n"
	print "//////////////////////////////////////////"
	print "Finished checking runningCouplings module!"
	print "//////////////////////////////////////////"