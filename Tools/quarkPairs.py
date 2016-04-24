####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module for creating quark pairs to use during dipole showering."""

##Import required modules:##
import math
import random
import assertions
import precision
import constants
import particleData
import runningCouplings
import kinematics

print "\n//////////////////////////"
print "Loading quarkPairs module:"
print "Test code not written!"
print "//////////////////////////\n"
assertions.pause(__name__) #Can remove import above lower when removed.

##Initiate coupling constants:##
oneLoopAlphaEM = runningCouplings.oneLoopAlphaEM()

##Functions:##

def calc_chi(S123):
	"""A function to return a value of the function chi used for fermion tree-level differential cross-sections."""
	##Using p443 of 'Dynamics of the Standard Model' textbook.
	assert (type(S123) == float)
	__zCode = particleData.knownParticles.get_code_from_name('Z-boson')
	__mZ = particleData.knownParticles.get_mass_from_code(__zCode)
	__gammaZ = particleData.knownParticles.get_width_from_code(__zCode)
	__sWSquared = constants.sin_squared_thetaW()
	__sW, __cW = math.sqrt(__sWSquared), math.sqrt(1.0 - __sWSquared)
	__term1 = (1.0/(2.0*__sW*__cW))*(1.0/(2.0*__sW*__cW))
	__term2 = S123/(S123 - (__mZ*__mZ) + ((1.0j)*__mZ*__gammaZ))
	return __term1*__term2

def calc_mod_squared_chi(S123):
	"""A function to return mod squared of the function chi used for fermion tree-level differential cross-sections."""
	##Using p443 of 'Dynamics of the Standard Model' textbook.
	assert (type(S123) == float)
	__chi = calc_chi(S123)
	return (__chi*__chi.conjugate()).real ##The .real just removes the 0j.

def calc_real_chi(S123):
	"""A function to return mod squared of the function chi used for fermion tree-level differential cross-sections."""
	##Using p443 of 'Dynamics of the Standard Model' textbook.
	assert (type(S123) == float)
	__chi = calc_chi(S123)
	return __chi.real

def calc_f_for_cos_theta(S123,fermionCode,cosTheta):
	"""A function to return the tree-level cross-section w.r.t cos(theta) for a given fermion f in e+e- -> ffBar."""
	##Using p443 of 'Dynamics of the Standard Model' textbook.
	##Already integrated dphi as no phi dependence.
	assert (type(S123) == float)
	assert (type(fermionCode) == int)
	assert (particleData.knownParticles.is_fermion(fermionCode))
	assert ((type(cosTheta) == float) and (-1 <= cosTheta) and (cosTheta <= 1))
	__alphaEMNow = oneLoopAlphaEM.calculate(S123)
	__electronCode = particleData.knownParticles.get_code_from_name('electron')
	__cNWC = constants.neutral_weak_couplings
	__gVE, __gAE = __cNWC('v',__electronCode), __cNWC('a',__electronCode)
	__gVF, __gAF = __cNWC('v',fermionCode), __cNWC('a',fermionCode)
	__mainPrefactor = (math.pi*__alphaEMNow*__alphaEMNow)/(2.0*S123)
	__term1A = (1.0 + (2.0*__gVE*__gVF*calc_real_chi(S123)))	
	__term1B = ((__gVE*__gVE) + (__gAE*__gAE))*((__gVF*__gVF) + (__gAF*__gAF))*calc_mod_squared_chi(S123)
	__term1C = (1.0 + (cosTheta*cosTheta))
	__term1 = (__term1A + __term1B)*__term1C
	__term2A = 4.0*__gAE*__gAF*calc_real_chi(S123)
	__term2B = 8.0*__gVE*__gAE*__gVF*__gAF*calc_mod_squared_chi(S123)
	__term2 = (__term2A + __term2B)*cosTheta
	__result = __mainPrefactor*(__term1 + __term2)
	return __result

def calc_g_for_cos_theta(S123,fermionCode,cosTheta):
	"""A function to return the primitive function of the overestimation in the theta sudakov."""
	##Using p443 of 'Dynamics of the Standard Model' textbook.
	##Already integrated dphi as no phi dependence.
	assert (type(S123) == float)
	assert (type(fermionCode) == int)
	assert (particleData.knownParticles.is_fermion(fermionCode))
	assert ((type(cosTheta) == float) and (-1 <= cosTheta) and (cosTheta <= 1))
	__alphaEMNow = oneLoopAlphaEM.calculate(S123)
	__electronCode = particleData.knownParticles.get_code_from_name('electron')
	__cNWC = constants.neutral_weak_couplings
	__gVE, __gAE = __cNWC('v',__electronCode), __cNWC('a',__electronCode)
	__gVF, __gAF = __cNWC('v',fermionCode), __cNWC('a',fermionCode)
	__mainPrefactor = (math.pi*__alphaEMNow*__alphaEMNow)/(2.0*S123)
	__term1A = (1.0 + (2.0*__gVE*__gVF*calc_real_chi(S123)))	
	__term1B = ((__gVE*__gVE) + (__gAE*__gAE))*((__gVF*__gVF) + (__gAF*__gAF))*calc_mod_squared_chi(S123)
	##Over-estimate with 1.2*__term1 for term1+term2 (inc cos fns.) and then 2 for (1.0 + cosx^2).
	__result = 1.2*2.0*__mainPrefactor*(__term1A + __term1B)
	return __result

def calc_G_for_cos_theta(S123,fermionCode,cosTheta):
	"""A function to return the primitive function of the overestimation in the theta sudakov."""
	##Using p443 of 'Dynamics of the Standard Model' textbook.
	##Already integrated dphi as no phi dependence.
	assert (type(S123) == float)
	assert (type(fermionCode) == int)
	assert (particleData.knownParticles.is_fermion(fermionCode))
	assert ((type(cosTheta) == float) and (-1 <= cosTheta) and (cosTheta <= 1))
	__alphaEMNow = oneLoopAlphaEM.calculate(S123)
	__electronCode = particleData.knownParticles.get_code_from_name('electron')
	__cNWC = constants.neutral_weak_couplings
	__gVE, __gAE = __cNWC('v',__electronCode), __cNWC('a',__electronCode)
	__gVF, __gAF = __cNWC('v',fermionCode), __cNWC('a',fermionCode)
	__mainPrefactor = (math.pi*__alphaEMNow*__alphaEMNow)/(2.0*S123)
	__term1A = (1.0 + (2.0*__gVE*__gVF*calc_real_chi(S123)))	
	__term1B = ((__gVE*__gVE) + (__gAE*__gAE))*((__gVF*__gVF) + (__gAF*__gAF))*calc_mod_squared_chi(S123)
	##Over-estimate with 1.2*__term1 for term1+term2 (inc cos fns.) and then 2 for (1.0 + cosx^2).
	__result = 1.2*2.0*__mainPrefactor*(__term1A + __term1B)*cosTheta
	return __result

def calc_inverse_G_for_cos_theta(S123,fermionCode,valueIn):
	"""A function to return the inverse of the primitive function of the overestimation in the theta sudakov."""
	##Using p443 of 'Dynamics of the Standard Model' textbook.
	##Already integrated dphi as no phi dependence.
	assert (type(S123) == float)
	assert (type(fermionCode) == int)
	assert (particleData.knownParticles.is_fermion(fermionCode))
	assert (type(valueIn) == float)
	__alphaEMNow = oneLoopAlphaEM.calculate(S123)
	__electronCode = particleData.knownParticles.get_code_from_name('electron')
	__cNWC = constants.neutral_weak_couplings
	__gVE, __gAE = __cNWC('v',__electronCode), __cNWC('a',__electronCode)
	__gVF, __gAF = __cNWC('v',fermionCode), __cNWC('a',fermionCode)
	__mainPrefactor = (math.pi*__alphaEMNow*__alphaEMNow)/(2.0*S123)
	__term1A = (1.0 + (2.0*__gVE*__gVF*calc_real_chi(S123)))	
	__term1B = ((__gVE*__gVE) + (__gAE*__gAE))*((__gVF*__gVF) + (__gAF*__gAF))*calc_mod_squared_chi(S123)
	##Over-estimate with 1.2*__term1 for term1+term2 (inc cos fns.) and then 2 for (1.0 + cosx^2).
	__result = valueIn/(1.2*2.0*__mainPrefactor*(__term1A + __term1B))
	return __result

def get_quarks_theta(S123,fermionCode):
	"""A function to get theta for a q-qBar pair produced from an e+e- collision."""
	##Using veto algorithm in "PYTHIA 6.0 Physics and Manual".
	assert (type(S123) == float)
	assert (type(fermionCode) == int)
	assert (fermionCode in particleData.knownParticles.get_known_quarks())
	__cosThetaMin, __cosThetaMax = -1.0, 1.0
	__notChosen1 = True
	while __notChosen1:
		__R1 = random.random()
		__Gmin = calc_G_for_cos_theta(S123,fermionCode,__cosThetaMin)
		__Gmax = calc_G_for_cos_theta(S123,fermionCode,__cosThetaMax)
		__cosThetai = calc_inverse_G_for_cos_theta(S123,fermionCode,(__R1*(__Gmax - __Gmin) + __Gmin))
		assert ((-1.0 < __cosThetai) and (__cosThetai < 1.0))
		__R2 = random.random()
		__ratio = calc_f_for_cos_theta(S123,fermionCode,__cosThetai) / calc_g_for_cos_theta(S123,fermionCode,__cosThetai)
		if (__ratio <= __R2):
			continue ##Don't accept, notChosen1 == True, while will run again.
		else:
			return math.acos(__cosThetai)

def get_quark_weights(S123,posQuarkCodes,cosTheta = None):
	"""A function to calculate the quark weights to beused by the ME generator."""
	__cE = constants.cut_off_energy()
	assert ((type(S123) == float) and (S123 > __cE*__cE))
	assert (type(posQuarkCodes) == list)
	assert (((type(cosTheta) == float) and (-1.0 <= cosTheta) and (cosTheta <= 1.0)) or (cosTheta == None))
	##Using dictionary not list so you don't have to assume an order of the quark codes.
	__hardCodedQWeights = {1:16.05764282,2:11.94029851,3:16.05764282,4:24.32321153,5:31.62120432}
	__qWeights = {} ##Empty dictionary initiated.
	for __qCode in posQuarkCodes: ##Uses the quark codes given to weight for the chosen set.
		if (cosTheta == None):
			__qWeights[__qCode] = __hardCodedQWeights[__qCode] ##Only take the hard coded ones required so normalisation works.
		else:
			__qWeights[__qCode] = calc_f_for_cos_theta(S123,__qCode,cosTheta)
	__total = sum(__qWeights.itervalues())
	__qWeights.update((x, y/__total) for x, y in __qWeights.items()) ##Normalise them
	assert precision.check_numbers_equal(sum(__qWeights.itervalues()),1.0) ##i.e they should now be normalised to 1.
	return __qWeights

def get_quark_code(S123,posQuarkCodes,cosTheta = None):
	"""A function to return a weighted random quark code."""
	__cE = constants.cut_off_energy()
	assert ((type(S123) == float) and (S123 > __cE*__cE))
	assert (type(posQuarkCodes) == list)
	assert (((type(cosTheta) == float) and (-1.0 <= cosTheta) and (cosTheta <= 1.0)) or (cosTheta == None))
	__qWeights = get_quark_weights(S123,posQuarkCodes,cosTheta) ##Re-set for every cosTheta.
	__R = random.random() ##Doesn't include 1 so using < not <= for checks is fair.
	__sumSoFar = 0.0
	for __qCode in posQuarkCodes: ##Using the initialised set for which the weights are normalised to 1.
		__sumSoFar += __qWeights[__qCode]
		if (__R < __sumSoFar):
			return __qCode
	assert False ##Should have returned by now!

def get_code_and_theta(S123,posQuarkCodes):
	"""A function to select a quark code and cos(theta) value consistent with the differential cross section."""
	__cE = constants.cut_off_energy()
	assert ((type(S123) == float) and (S123 > __cE*__cE))
	assert (type(posQuarkCodes) == list)
	__qCode = get_quark_code(S123,posQuarkCodes)
	__theta = get_quarks_theta(S123,__qCode)
	return __qCode, __theta

##Classes:##

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##
	import matplotlib.ticker as mtick
	from matplotlib import pyplot
	import numpy

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "//////////////////////////"
	print "Testing quarkPairs module:"
	print "//////////////////////////"
	assertions.pause(__name__)
	
	##Setup here:##
	zCode = particleData.knownParticles.get_code_from_name('Z-boson')
	Mz = particleData.knownParticles.get_mass_from_code(zCode)
	numberThetas = 100000
	thetaRange = [i*math.pi/(numberThetas - 1.0) for i in range(0,numberThetas)]
	cosThetaRange = [math.cos(i*math.pi/(numberThetas - 1.0)) for i in range(0,numberThetas)]
	testS123, testQCodes, testPosQCodes = Mz*Mz, [1,2,3,4,5,6], [1,2,3]
	numTestIts = float(numberThetas)
	numberPerBin = 100
	numberBins = int(numTestIts)/numberPerBin
	testQNames = [r"$\rm{d-quark}$",r"$\rm{u-quark}$",r"$\rm{s-quark}$",r"$\rm{c-quark}$",r"$\rm{b-quark}$",r"$\rm{t-quark}$"]
	numberS123s = 1000
	testS123s = [i*testS123/(numberS123s-1.0) for i in range(0,numberS123s)][1:] ##Remove 0 where undefined.

	###Test calc_real_chi:##
	#print "\n--------------------------------------------------\n"
	#print "Testing calc_real_chi:\n"
	#realChis = [calc_real_chi(anS123) for anS123 in testS123s]
	#pyplot.figure()
	#pyplot.title(r"$Re(\chi)\ as\ a\ function\ of\ S_{123}$")
	#pyplot.xlabel(r"$S_{123} (GeV^{2}/c^{4})$")
	#pyplot.ylabel(r"$Re(\chi(S_{123}))$")
	#pyplot.plot(testS123s,realChis,linewidth = 2, label = r"$Re(\chi(S_{123}))$")
	#pyplot.legend()
	#assertions.show_graph()
	#print "\nFinished testing calc_real_chi."
	#assertions.pause(__name__)
#
	###Test calc_mod_squared_chi:##
	#print "\n--------------------------------------------------\n"
	#print "Testing calc_mod_squared_chi:\n"
	#realChis = [calc_mod_squared_chi(anS123) for anS123 in testS123s]
	#pyplot.figure()
	#pyplot.title(r"$\|\chi\|^{2}\ as\ a\ function\ of\ S_{123}$")
	#pyplot.xlabel(r"$S_{123} (GeV^{2}/c^{4})$")
	#pyplot.ylabel(r"$\|\chi(S_{123})\|^{2}$")
	#pyplot.plot(testS123s,realChis,linewidth = 2, label = r"$\|\chi(S_{123})\|^{2}$")
	#pyplot.legend()
	#assertions.show_graph()
	#print "\nFinished testing calc_mod_squared_chi."
	#assertions.pause(__name__)
#
	###Test prefactors in calc_f_for_cos_theta:##
	#print "\n--------------------------------------------------\n"
	#print "Testing prefactors in calc_f_for_cos_theta:\n"
	#alphaEM = oneLoopAlphaEM.calculate
	#electronCode = particleData.knownParticles.get_code_from_name('electron')
	#cNWC = constants.neutral_weak_couplings
	#for fCode in testQCodes:
	#	gVE, gAE = cNWC('v',electronCode), cNWC('a',electronCode)
	#	gVF, gAF = cNWC('v',fCode), cNWC('a',fCode)
	#	prefactors = [(math.pi*alphaEM(anS123)*alphaEM(anS123))/(2.0*anS123) for anS123 in testS123s]
	#	term1As = [(1.0 + (2.0*gVE*gVF*calc_real_chi(anS123))) for anS123 in testS123s]
	#	term1Bs = [(((gVE*gVE) + (gAE*gAE))*((gVF*gVF) + (gAF*gAF))*calc_mod_squared_chi(anS123)) for anS123 in testS123s]
	#	term1s = [(term1As[i] + term1Bs[i]) for i,anS123 in enumerate(testS123s)]
	#	term2As = [(4.0*gAE*gAF*calc_real_chi(anS123)) for anS123 in testS123s]
	#	term2Bs = [(8.0*gVE*gAE*gVF*gAF*calc_mod_squared_chi(anS123)) for anS123 in testS123s]
	#	term2s = [(term2As[i] + term2Bs[i]) for i,anS123 in enumerate(testS123s)]
	#	greaterThan = True
	#	for i, x in enumerate(term1s):
	#		if term2s[i] >= 2*term1s[i]:
	#			greaterThan = False
	#	if (greaterThan == True):
	#		print "Using " + str(fCode) + "; 2*term1 > term2 for all!"
	#	else:
	#		print "g failed!"
	#	pyplot.figure()
	#	pyplot.title(r"$Prefactors(S_{123})\ in\ fermion\_cross\_section\_\rm{d}\theta\ for\ $" + testQNames[fCode-1] + r"$s$")
	#	pyplot.xlabel(r"$S_{123} (GeV^{2}/c^{4})$")
	#	pyplot.ylabel(r"$Prefactor(S_{123})$")
	#	pyplot.plot(testS123s,prefactors,linewidth = 2, label = r"$Main\ prefactor$")
	#	pyplot.plot(testS123s,term1s,linewidth = 2, label = r"$Prefactor\ 1$")
	#	pyplot.plot(testS123s,term2s,linewidth = 2, label = r"$Prefactor\ 2$")
	#	pyplot.legend()
	#	assertions.show_graph()
	#print "\nFinished testing prefactors in calc_f_for_cos_theta."
	#assertions.pause(__name__)
#
	###Test terms in calc_f_for_cos_theta:##
	#print "\n--------------------------------------------------\n"
	#print "Testing terms in calc_f_for_cos_theta:\n"
	#alphaEM = oneLoopAlphaEM.calculate
	#electronCode = particleData.knownParticles.get_code_from_name('electron')
	#cNWC = constants.neutral_weak_couplings
	#for fCode in testQCodes:
	#	gVE, gAE = cNWC('v',electronCode), cNWC('a',electronCode)
	#	gVF, gAF = cNWC('v',fCode), cNWC('a',fCode)
	#	prefactors = [(math.pi*alphaEM(testS123)*alphaEM(testS123))/(2.0*testS123) for aCosTheta in cosThetaRange]
	#	term1As = [(1.0 + (2.0*gVE*gVF*calc_real_chi(testS123))) for aCosTheta in cosThetaRange]
	#	term1Bs = [(((gVE*gVE) + (gAE*gAE))*((gVF*gVF) + (gAF*gAF))*calc_mod_squared_chi(testS123)) for aCosTheta in cosThetaRange]
	#	term1s = [(term1As[i] + term1Bs[i])*(1.0 + aCosTheta*aCosTheta) for i,aCosTheta in enumerate(cosThetaRange)]
	#	term2As = [(4.0*gAE*gAF*calc_real_chi(testS123)) for aCosTheta in cosThetaRange]
	#	term2Bs = [(8.0*gVE*gAE*gVF*gAF*calc_mod_squared_chi(testS123)) for aCosTheta in cosThetaRange]
	#	term2s = [(term2As[i] + term2Bs[i])*aCosTheta for i,aCosTheta in enumerate(cosThetaRange)]
	#	differences = [((2*term1s[i]) - term2s[i]) for i,aCosTheta in enumerate(cosThetaRange)]
	#	greaterThan = True
	#	for i, x in enumerate(term1s):
	#		if term1s[i] + term2s[i] >= 1.2*term1s[i]:
	#			greaterThan = False
	#	if (greaterThan == True):
	#		print "Using " + str(fCode) + "; 1.2*term1 > term1 + term2 for all!"
	#	else:
	#		print "g failed!"
	#	pyplot.figure()
	#	pyplot.title(r"$Terms\ in\ fermion\_cross\_section\_\rm{d}\theta\ for\ $" + testQNames[fCode-1] + r"$s$")
	#	pyplot.xlabel(r"$S_{123} (GeV^{2}/c^{4})$")
	#	pyplot.ylabel(r"$Prefactor(S_{123})$")
	#	pyplot.plot(cosThetaRange,prefactors,linewidth = 2, label = r"$Main\ prefactor$")
	#	pyplot.plot(cosThetaRange,term1s,linewidth = 2, label = r"$Prefactor\ 1$")
	#	pyplot.plot(cosThetaRange,term2s,linewidth = 2, label = r"$Prefactor\ 2$")
	#	pyplot.plot(cosThetaRange,differences,linewidth = 2, label = r"$2\times(1)\ -\ (2)$")
	#	pyplot.legend()
	#	assertions.show_graph()
	#print "\nFinished testing terms in calc_f_for_cos_theta."
	#assertions.pause(__name__)
#
	###Test calc_f_for_cos_theta and calc_g_for_cos_theta:##
	#print "\n--------------------------------------------------\n"
	#print "Testing calc_f_for_cos_theta and calc_g_for_cos_theta:\n"
	#cosThetafValues, cosThetagValues = [], []
	##combine into one graph?
	#for qCode in testQCodes:
	#	cosThetafValues.append([calc_f_for_cos_theta(testS123,qCode,i) for i in cosThetaRange])
	#	cosThetagValues.append([calc_g_for_cos_theta(testS123,qCode,i) for i in cosThetaRange])
	#	##Don't want to normalise here as just looking at raw f and g.
	#	greaterThan = True
	#	for i, x in enumerate(cosThetafValues):
	#		if cosThetafValues[i] >= cosThetagValues[i]:
	#			greaterThan = False
	#	if (greaterThan == True):
	#		print "Using " + qCode + "; g > f for all!"
	#	else:
	#		print "g failed!"
	#	pyplot.figure()
	#	pyplot.title(r"$Comparing\ the\ functions\ f\ and\ g\ for\ $" + testQNames[qCode - 1] + r"$s$")
	#	pyplot.xlabel(r"$\theta\ (rad)$")
	#	pyplot.ylabel(r"$Normalised\ \sigma(\theta)$")
	#	pyplot.ylim(0.0,0.000015)
	#	pyplot.xlim(-1.1,1.1)
	#	pyplot.axvline(-1, linestyle = '--', color = 'red')
	#	pyplot.axvline(1, linestyle = '--', color = 'red')
	#	pyplot.plot(cosThetaRange,cosThetafValues[qCode - 1],linewidth = 2, label = r"$f(cos\theta)$")
	#	pyplot.plot(cosThetaRange,cosThetagValues[qCode - 1],linewidth = 2, label = r"$g(cos\theta)$")
	#	pyplot.legend()
	#	assertions.show_graph()
	#print "\nFinished testing calc_f_for_cos_theta and calc_g_for_cos_theta."
	#assertions.pause(__name__)
#
	###Test get_quarks_theta:##
	#print "\n--------------------------------------------------\n"
	#print "Testing get_quarks_theta:\n"
	#cosThetaProbs = []
	#allCosThetaHists, allCosThetaBins, allBinCentres = [], [], []
	#cosThetasOut = []
	#for qCode in testQCodes:
	#	print "Using" + str(qCode[4:-1]) +":"
	#	currentCosThetaProbs = [calc_f_for_cos_theta(testS123,qCode,aCosTheta) for aCosTheta in cosThetaRange]
	#	thetaProbsSum = sum(currentCosThetaProbs)
	#	cosThetaProbs.append([i/float(thetaProbsSum) for i in currentCosThetaProbs])
	#	print "Normalised expected values sum to:", sum(cosThetaProbs[-1])
	#	cosThetasOut.append([math.cos(get_quarks_theta(testS123,qCode)) for n in range(int(numTestIts))])
	#	bins = numpy.linspace(-1.0,1.0,numberBins+1)
	#	cosThetasHists, cosThetaBins = numpy.histogram(cosThetasOut[-1],bins)
	#	##Normalise to 1. Have to divide by the number per bin for plotting the average.
	#	allCosThetaHists.append([aNumberInBin/float(numTestIts*numberPerBin) for aNumberInBin in cosThetasHists])
	#	print "Normalised generated values sum to:", sum(allCosThetaHists[-1])*numberPerBin ##As divided each by above.
	#	allCosThetaBins.append(cosThetaBins)
	#	allBinCentres.append([])
	#	for i in range(len(allCosThetaBins[-1])-1): ##-1 for number of centres
	#		allBinCentres[-1].append((allCosThetaBins[-1][i] + allCosThetaBins[-1][i+1])/2.0)
	#xaxes = [r"$\rm{cos\ \theta}$",r"$\rm{cos\ \theta}$",r"$\rm{cos\ \theta}$",r"$\rm{cos\ \theta}$",r"$\rm{cos\ \theta}$",r"$\rm{cos\ \theta}$"]
	#yaxes = [r"$\rm{P}(\rm{cos\ \theta})$",r"$\rm{P}(\rm{cos\ \theta})$",r"$\rm{P}(\rm{cos\ \theta})$",r"$\rm{P}(\rm{cos\ \theta})$"]
	#yaxes += [r"$\rm{P}(\rm{cos\ \theta})$",r"$\rm{P}(\rm{cos\ \theta})$"]
	#yMaxs = [0.0025,0.0025,0.0025,0.0025,0.0025,0.0025]
	#supTitle = r"$\rm{Monte\ Carlo\ sampling\ of\ cos\ \theta\ for}\ e^{+}e^{-}\ \rightarrow\ q\bar{q}\ \rm{using\ "
	#supTitle += str(int(numTestIts)) + r"\ iterations}$"
	#figure,axes = pyplot.subplots(3,2)
	#axes = axes.ravel()
	#for idx,ax in enumerate(axes):
	#	ax.plot(allBinCentres[idx], allCosThetaHists[idx], linestyle = "solid", color = "blue", linewidth = 2)
	#	ax.plot(cosThetaRange, cosThetaProbs[idx], linestyle = "solid", color = "red", linewidth = 2)
	#	ax.set_title(testQNames[idx])
	#	ax.set_xlabel(xaxes[idx])
	#	ax.set_ylabel(yaxes[idx])
	#	ax.axis([-1,1,0,yMaxs[idx]])
	#	ax.set_yticks([0.0,0.001,0.002,0.003],minor=False)
	#	ax.yaxis.set_major_formatter(mtick.FixedFormatter([r"$0.0$",r"$1.0$",r"$2.0$",r"$\times\ 10^{-3}$"]))
	#	ax.set_xticks([-1.0,-0.5,0.0,0.5,1.0],minor=False)
	#	ax.xaxis.set_major_formatter(mtick.FixedFormatter([r"$\minus1.0$",r"$\minus0.5$",r"$0.0$",r"$0.5$",r"$1.0$"]))
	#pyplot.suptitle(supTitle,fontsize = "16")
	#line1 = pyplot.Line2D((0,1),(0,0), color="blue", linewidth = 2)
	#line2 = pyplot.Line2D((0,1),(0,0), color="red", linewidth = 2)
	#lines, figStrs = [line1,line2], [r"$\rm{Generated}$",r"$\rm{Expected}$"]
	#figure.legend(lines, figStrs, bbox_to_anchor=[0.5, 0.05],loc='center', ncol=2)
	#pyplot.tight_layout()
	###Space main title out to prevent overlapping and allow space for legend below:
	#pyplot.subplots_adjust(top=0.85,bottom=0.17)
	#assertions.show_graph()
	#print "\nFinished testing get_quarks_theta."
	#assertions.pause(__name__)

	##Test get_quark_weights and get_quark_code:##
	print "\n--------------------------------------------------\n"
	print "Testing get_quark_weights and get_quark_code:\n"
	testQuarkCodes = [get_quark_code(testS123,testPosQCodes) for i in range(int(numTestIts))]
	numsResults = [sum([i for i in testQuarkCodes if i == qCode]) for qCode in [testPosQCodes]]
	print len(testQuarkCodes), sum(numsResults)
	print numsResults
	sumPercent = 0.0
	for x, qCode in enumerate(testPosQCodes):
		sumPercent += numsResults[x]*100/numTestIts
		print "For code", qCode, "produced" ,numsResults[x]*100/numTestIts, "%"
	print "These add up to", sumPercent, "%"
	assertions.show_graph()
	print "\nFinished testing get_quark_weights and get_quark_code."
	assertions.pause(__name__)

	##Done testing:##
	print "\n---------------------------------------------\n"
	print "////////////////////////////////////"
	print "Finished checking quarkPairs module!"
	print "////////////////////////////////////"