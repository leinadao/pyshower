####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module for handling dipole showering gluon emission Sudakov form factors."""

##Using cross sections from the Ariadne paper throughout.##

##Import required modules:##
import random
import math
import assertions
import precision
import constants
import particleData
import runningCouplings
import kinematics

print "\n///////////////////////"
print "Loading sudakov module:"
print "Test code not re-written!"
print "///////////////////////\n"
assertions.pause_loading_module()

##Initiate coupling constants:##
alphaS = runningCouplings.oneLoopAlphaS()
alphaEM = runningCouplings.oneLoopAlphaEM()

##Functions:##

def check_rapidity_allowed(S123,PperpSquared,y):
	"""A function to check if the rapidity is in the allowed region of phase space."""
	##Using equation 49 of Initial-state showering based on colour dipoles connected to incoming parton lines.
	assert (type(S123) == float)
	assert (type(PperpSquared) == float)
	assert (type(y) == float)
	__yMagMax = math.acosh(math.sqrt(S123/PperpSquared)/2.0)
	if (abs(y) > __yMagMax):
		return False
	else:
		return True

##~~Cross Sections~~##

def qqBar_to_qgqBar_cS(S123,PperpSquared,y): ##difCrossSecId = 0
	"""A function to calculate the differential cross section for qqBar -> qgqBar."""
	assert assertions.all_are_numbers([S123,PperpSquared,y])
	__alphaSNow = alphaS.calculate(PperpSquared)
	__x1 = kinematics.calculate_x1(S123,math.sqrt(PperpSquared),y)
	__x3 = kinematics.calculate_x3(S123,math.sqrt(PperpSquared),y)
	__preFactor = (2.0*__alphaSNow)/(3.0*math.pi)
	__numerator = (__x1*__x1) + (__x3*__x3)
	__denominator = (1.0 - __x1)*(1.0 - __x3)
	return __preFactor*(__numerator/__denominator)

def qg_to_qgg_cS(S123,PperpSquared,y): ##difCrossSecId = 1
	"""A function to calculate the differential cross section for qg -> qgg."""
	assert assertions.all_are_numbers([S123,PperpSquared,y])
	__alphaSNow = alphaS.calculate(PperpSquared)
	__x1 = kinematics.calculate_x1(S123,math.sqrt(PperpSquared),y)
	__x3 = kinematics.calculate_x3(S123,math.sqrt(PperpSquared),y)
	__preFactor = (3.0*__alphaSNow)/(4.0*math.pi)
	__numerator = (__x1*__x1) + (__x3*__x3*__x3)
	__denominator = (1.0 - __x1)*(1.0 - __x3)
	return __preFactor*(__numerator/__denominator)

def gg_to_ggg_cS(S123,PperpSquared,y): ##difCrossSecId = 2
	"""A function to calculate the differential cross section for gg -> ggg."""
	assert assertions.all_are_numbers([S123,PperpSquared,y])
	__alphaSNow = alphaS.calculate(PperpSquared)
	__x1 = kinematics.calculate_x1(S123,math.sqrt(PperpSquared),y)
	__x3 = kinematics.calculate_x3(S123,math.sqrt(PperpSquared),y)
	__preFactor = (3.0*__alphaSNow)/(4.0*math.pi)
	__numerator = (__x1*__x1*__x1) + (__x3*__x3*__x3)
	__denominator = (1.0 - __x1)*(1.0 - __x3)
	return __preFactor*(__numerator/__denominator)

def qg_to_qQQBar_cS(S123,PperpSquared,y): ##difCrossSecId = 3
	"""A function to calculate the differential cross section for qg -> qQQBar."""
	assert assertions.all_are_numbers([S123,PperpSquared,y])
	__alphaSNow = alphaS.calculate(PperpSquared)
	__x1 = kinematics.calculate_x1(S123,math.sqrt(PperpSquared),y)
	__x2 = kinematics.calculate_x2(S123,math.sqrt(PperpSquared),y)
	__x3 = kinematics.calculate_x3(S123,math.sqrt(PperpSquared),y)
	__preFactor = (3.0*__alphaSNow)/(8.0*math.pi)
	__numerator = ((1.0 - __x1)*(1.0 - __x1)) + ((1.0 - __x2)*(1.0 - __x2))
	__denominator = (1.0 - __x3)
	return __preFactor*(__numerator/__denominator)

def gg_to_gQQBar_cS(S123,PperpSquared,y): ##difCrossSecId = 4,5
	"""A function to calculate the differential cross section for gg -> gQQBar."""
	assert assertions.all_are_numbers([S123,PperpSquared,y])
	__alphaSNow = alphaS.calculate(PperpSquared)
	__x1 = kinematics.calculate_x1(S123,math.sqrt(PperpSquared),y)
	__x2 = kinematics.calculate_x2(S123,math.sqrt(PperpSquared),y)
	__x3 = kinematics.calculate_x3(S123,math.sqrt(PperpSquared),y)
	__preFactor = (3.0*__alphaSNow)/(8.0*math.pi)
	__numerator = ((1.0 - __x1)*(1.0 - __x1)) + ((1.0 - __x2)*(1.0 - __x2))
	__denominator = (1.0 - __x3)
	return __preFactor*(__numerator/__denominator)

def qqBar_to_qpqBar_cS(S123,PperpSquared,y,code1,code2): ##difCrossSecId = 6
	"""A function to calculate the differential cross section for qqBar -> qpqBar."""
	assert assertions.all_are_numbers([S123,PperpSquared,y])
	__kQs = particleData.knownParticles.get_known_quarks()
	__expectedCodes = __kQs + [-x for x in __kQs]
	assert ((code1 in __expectedCodes) and (code2 in __expectedCodes))
	__alphaEMNow = alphaEM.calculate(PperpSquared)
	__charge1 = particleData.knownParticles.get_charge_from_code(code1)
	__charge2 = particleData.knownParticles.get_charge_from_code(code2)
	__x1 = kinematics.calculate_x1(S123,math.sqrt(PperpSquared),y)
	__x3 = kinematics.calculate_x3(S123,math.sqrt(PperpSquared),y)
	__preFactor = (__alphaEMNow*abs(__charge1)*abs(__charge2))/(2.0*math.pi)
	__numerator = (__x1*__x1) + (__x3*__x3)
	__denominator = (1.0 - __x1)*(1.0 - __x3)
	return __preFactor*(__numerator/__denominator)

def sum_cSs(difCrossSecIds,S123,PperpSquared,y,code1,code2):
	"""A function to return the sum of each cross section from its process code in the list given."""
	assert (type(difCrossSecIds) == list)
	assert assertions.all_are_numbers([S123,PperpSquared,y])
	__kQs = particleData.knownParticles.get_known_quarks()
	__gC = particleData.knownParticles.get_code_from_name('gluon')
	__expectedCodes = __kQs + [-x for x in __kQs] + [__gC]
	assert ((code1 in __expectedCodes) and (code2 in __expectedCodes))
	__result = 0.0
	for __difCrossSecId in difCrossSecIds:
		assert ((0 <= __difCrossSecId) and (__difCrossSecId <= 6))
		if __difCrossSecId == 0:
			__result += qqBar_to_qgqBar_cS(S123,PperpSquared,y)
		elif __difCrossSecId == 1:
			__result += qg_to_qgg_cS(S123,PperpSquared,y)
		elif __difCrossSecId == 2:
			__result += gg_to_ggg_cS(S123,PperpSquared,y)
		elif __difCrossSecId == 3:
			__result += qg_to_qQQBar_cS(S123,PperpSquared,y)
		elif __difCrossSecId == 4:
			__result += gg_to_gQQBar_cS(S123,PperpSquared,y)
		elif __difCrossSecId == 5:
			__result += gg_to_gQQBar_cS(S123,PperpSquared,y)
		elif __difCrossSecId == 6:
			__result += qqBar_to_qpqBar_cS(S123,PperpSquared,y,code1,code2)
	return __result

##~~Overestimate Functions~~##

def calc_g(PperpSquared,y): ##Calculate the function g(Pperpsquared,y) = g1(Pperpsquared)*g2(y).
	"""A function to calculate the overestimation function g for all processes."""
	assert assertions.all_are_numbers([PperpSquared,y])
	__alphaSMax = alphaS.get_shower_max()
	return (3.0*__alphaSMax) / (2.0*math.pi*PperpSquared)

def calc_G(PperpSquared,S123):
	"""A function to calculate the overestimation function G, the primitive of g(PperpSquared) = g1(PperpSquared)*int{g2(y) dy}, for all processes."""
	assert assertions.all_are_numbers([PperpSquared,S123])
	__alphaSMax = alphaS.get_shower_max()
	__term1 = ((-3.0)*__alphaSMax) / (4.0*math.pi)
	return __term1*math.log(S123/PperpSquared)*math.log(S123/PperpSquared)

def calc_inv_G(valueIn,S123):
	"""A function to calculate the inverse of G for for all processes."""
	assert assertions.all_are_numbers([valueIn,S123])
	__alphaSMax = alphaS.get_shower_max()
	__term1 = ((-4.0)*math.pi*valueIn) / (3.0*__alphaSMax)
	__exp = math.exp(math.sqrt(__term1))
	return S123/__exp

def calc_G2(y):
	"""A function to calculate the overestimation function G2, the primitive of g2, for all processes."""
	assert assertions.all_are_numbers([y])
	return y

def calc_inv_G2(valueIn):
	"""A function to calculate the inverse of G2 for for all processes."""
	assert assertions.all_are_numbers([valueIn])
	return valueIn

##~~Top Functions~~##

def solve_sudakovs(PperpSquaredMax,S123,difCrossSecIds,code1,code2):
	"""A function to generate a PperpSquared and y for an emission or splitting, using the Sudakov form factor."""
	##Codes 1,2 still needed to get the charges for EM radiation.
	##Using veto algorithm in "PYTHIA 6.0 Physics and Manual".
	##And p198 of "Monte Carlo simulations of hard QCD radiation" for bivariant algorithm.
	##And p16 of "Initial-state showering based on colour dipoles connected to incoming parton lines" for real y limits.
	##And p5 of "Fooling around with the Sudakov veto algorithm" for multiple emission form.
	assert assertions.all_are_numbers([PperpSquaredMax,S123])
	assert (type(difCrossSecIds) == list)
	__kQs = particleData.knownParticles.get_known_quarks()
	__gC = particleData.knownParticles.get_code_from_name('gluon')
	__expectedCodes = __kQs + [-x for x in __kQs] + [__gC]
	assert ((code1 in __expectedCodes) and (code2 in __expectedCodes))
	__cutOff = constants.cut_off_energy()*constants.cut_off_energy()
	if (PperpSquaredMax < __cutOff):
		return None, None, None
	__notChosen1, __y, __PperpSquaredi = True, 0, assertions.force_float_number(PperpSquaredMax)
	__PperpSquarediMinus1 = __PperpSquaredi
	while __notChosen1: ##i += 1
		__notChosen2 = True
		while __notChosen2:
			__R1 = random.random()
			__PperpSquaredi = calc_inv_G((math.log(__R1) + calc_G(__PperpSquarediMinus1,S123)),S123)
			if (__PperpSquaredi > S123/4.0): ##The maximum physically possible!
				__notChosen2 = True
				__PperpSquarediMinus1 = __PperpSquaredi
				continue #Re-run the while loop
			if (__PperpSquaredi <= __PperpSquarediMinus1): ##Don't update __PperpSquaredi if fails as would start allowing higher values!
				__notChosen2 = False ##Move on
		if (__PperpSquaredi < __cutOff):
			return None, None, None
		__R2 = random.random()
		__ymin, __ymax = -0.5*math.log(S123/__PperpSquaredi), 0.5*math.log(S123/__PperpSquaredi)
		__yi = calc_inv_G2(__R2*(calc_G2(__ymax) - calc_G2(__ymin)) + calc_G2(__ymin))
		if not check_rapidity_allowed(S123,__PperpSquaredi,__yi): ##Veto incorrect y-range. On the boundary is allowed.
			__PperpSquarediMinus1 = __PperpSquaredi
			continue ##Don't accept, notChosen1 == True, while will run again.
		__R3 = random.random()
		__numberCSIds = len(difCrossSecIds)
		__ratio = sum_cSs(difCrossSecIds,S123,__PperpSquaredi,__yi,code1,code2) / __numberCSIds*calc_g(__PperpSquaredi,__yi)
		if (__ratio <= __R3):
			__PperpSquarediMinus1 = __PperpSquaredi
			continue ##Don't accept, notChosen1 == True, while will run again.
		else:
			__notChosen1 = False ##Accept.
	##Prepare weights for choosing Id of process occuring:
	__idWeights = {} ##Empty dictionary initiated.
	for __difCrossSecId in difCrossSecIds: ##Uses the Ids given so will weight for the chosen set.
		__idWeights[__difCrossSecId] = sum_cSs([__difCrossSecId],S123,__PperpSquaredi,__yi,code1,code2)
	__total = sum(__idWeights.itervalues())
	__idWeights.update((x, y/__total) for x, y in __idWeights.items()) ##Normalise them
	assert precision.check_numbers_equal(sum(__idWeights.itervalues()),1.0) ##i.e they should now be normalised to 1.
	##Choose which process:
	__R4 = random.random() ##Doesn't include 1 -> using < not <= for checks is fair.
	__sumSoFar = 0.0
	__chosenId = None
	__chosen = False
	for __difCrossSecId in difCrossSecIds:
		if (not __chosen):
			__sumSoFar += __idWeights[__difCrossSecId]
			if (__R4 < __sumSoFar):
				__chosenId = __difCrossSecId
				__chosen = True
	assert (not (__chosenId == None)) ##Should have chosen by now!
	return __PperpSquaredi, __yi, __chosenId

def solve_case_group_0(PperpSquaredMax,S123,code1,code2):
	"""A function to return the next Pperp^2, y and winning process for case group 0/A."""
	##qqBar, qBarq, qq or qBarqBar so gluon emission or photon emission.
	assert assertions.all_are_numbers([PperpSquaredMax,S123])
	__kQs = particleData.knownParticles.get_known_quarks()
	__expectedCodes = __kQs + [-x for x in __kQs]
	assert ((code1 in __expectedCodes) and (code2 in __expectedCodes))
	__difCrossSecIds = [0] ##Gluon emission always possible.
	if assertions.EM_radiation_on():
		__difCrossSecIds.append(6) ##Photon emission possible.
	__PperpSquared, __y, __difCrossSecId = solve_sudakovs(PperpSquaredMax,S123,__difCrossSecIds,code1,code2)
	__convertIdToProcessCode = {None:0,0:1,1:1,2:1,6:3}
	return __PperpSquared, __y, __convertIdToProcessCode[__difCrossSecId]

def solve_case_group_1(PperpSquaredMax,S123,code1,code2):
	"""A function to return the next Pperp^2, y and winning process for case group 1/B."""
	##qqBar or qBarq so gluon emission or gluon splitting.
	assert assertions.all_are_numbers([PperpSquaredMax,S123])
	__kQs = particleData.knownParticles.get_known_quarks()
	__gC = particleData.knownParticles.get_code_from_name('gluon')
	__expectedCodes = __kQs + [-x for x in __kQs] + [__gC]
	assert ((code1 in __expectedCodes) and (code2 in __expectedCodes))
	assert ((code1 == __gC) or (code2 == __gC)) ##One must be a gluon.
	assert (not ((code1 == __gC) and (code2 == __gC))) ##One must be q/qBar.
	__difCrossSecIds = [1] ##Gluon emission always possible.
	if assertions.gluon_splitting_on():
		__difCrossSecIds.append(3) ##Gluon splitting possible.
	__PperpSquared, __y, __difCrossSecId = solve_sudakovs(PperpSquaredMax,S123,__difCrossSecIds,code1,code2)
	__convertIdToProcessCode = {None:0,0:1,1:1,2:1,3:2,4:2,5:2}
	return __PperpSquared, __y, __convertIdToProcessCode[__difCrossSecId]

def solve_case_group_2(PperpSquaredMax,S123,code1,code2):
	"""A function to return the next Pperp^2, y and winning process for case group 2/C."""
	##qqBar or qBarq so gluon emission, gluon 1 splitting or gluon 2 splitting.
	assert assertions.all_are_numbers([PperpSquaredMax,S123])
	__gC = particleData.knownParticles.get_code_from_name('gluon')
	assert ((code1 == __gC) and (code2 == __gC))
	__difCrossSecIds = [2] ##Gluon emission always possible.
	if assertions.gluon_splitting_on():
		__difCrossSecIds.append(4) ##Gluon 1 splitting possible.
		__difCrossSecIds.append(5) ##Gluon 2 splitting possible.
		##Note that gluons 1,2 to not refer to a specific gluon in the dipole, only that there are two gluons.
	__PperpSquared, __y, __difCrossSecId = solve_sudakovs(PperpSquaredMax,S123,__difCrossSecIds,code1,code2)
	__convertIdToProcessCode = {None:0,0:1,1:1,2:1,3:2,4:2,5:2}
	return __PperpSquared, __y, __convertIdToProcessCode[__difCrossSecId]

def solve(PperpSquaredMax,S123,code1,code2):
	"""A function to return the next Pperp^2, y and winning process for a given dipole."""
	assert assertions.all_are_numbers([PperpSquaredMax,S123])
	__kQs = particleData.knownParticles.get_known_quarks()
	__gC = particleData.knownParticles.get_code_from_name('gluon')
	__expectedCodes = __kQs + [-x for x in __kQs] + [__gC]
	assert ((code1 in __expectedCodes) and (code2 in __expectedCodes))
	if ((code1 in __kQs) and ((-code2) in __kQs)): ##qqBar.
		__caseGroupId = 0
	elif (((-code1) in __kQs) and (code2 in __kQs)): ##qBarq.
		__caseGroupId = 0
	elif (((code1) in __kQs) and (code2 in __kQs)): ##qq.
		__caseGroupId = 0
	elif (((-code1) in __kQs) and (-code2 in __kQs)): ##qBarqBar.
		__caseGroupId = 0
	elif ((code1 in __kQs) and (code2 == __gC)): ##qg.
		__caseGroupId = 1
	elif (((-code1) in __kQs) and (code2 == __gC)): ##qBarg.
		__caseGroupId = 1
	elif ((code1 == __gC) and (code2 in __kQs)): ##gq.
		__caseGroupId = 1
	elif ((code1 == __gC) and ((-code2) in __kQs)): ##gqBar.
		__caseGroupId = 1
	elif ((code1 == __gC) and (code2 == __gC)): ##gg.
		__caseGroupId = 2
	if (__caseGroupId == 0):
		__newPperpSquared, __newY, __newProcessCode = solve_case_group_0(PperpSquaredMax,S123,code1,code2)
	elif (__caseGroupId == 1):
		__newPperpSquared, __newY, __newProcessCode = solve_case_group_1(PperpSquaredMax,S123,code1,code2)
	elif (__caseGroupId == 2):
		__newPperpSquared, __newY, __newProcessCode = solve_case_group_2(PperpSquaredMax,S123,code1,code2)
 	##Process codes: 0 = stop shower, 1 = gluon emission, 2 = gluon splitting, 3 = photon emission.
 	##gg has two chances for a gluon to split but nothing specifies which is which so if occurs, randomly chosen later.
	return __newPperpSquared, __newY, __newProcessCode

##Classes:##

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##
	import numpy
	import matplotlib.pyplot as pyplot

#	##Define test functions:##
#	def calc_test_qqBar_int_func(PperpSquared,S123):
#		"""A function to return a value for the function in the Pperp^2 test integral, for qqBar."""
#		assert assertions.all_are_numbers([PperpSquared,S123])
#		alphaSNow = alphaS.calculate(PperpSquared)
#		term1 = (2.0*alphaSNow)/(3.0*math.pi*PperpSquared)
#		term2 = (2.0*math.log(S123/PperpSquared)) - 3.0 + ((4.0*PperpSquared)/S123) - ((PperpSquared*PperpSquared)/(S123*S123))
#		return term1*term2
#
#	def calc_test_qqBar_probs(S123,numtestIts):
#		"""A function to calculate the test probabilities as a function of Pperp^2, for qqBar."""
#		assert assertions.all_are_numbers([S123,numtestIts])
#		intCutOff = constants.cut_off_energy()*constants.cut_off_energy()
#		topLimit = S123 ##PperpSquaredMax = S123 for first dipole emission.
#		testValues = numpy.linspace(topLimit,intCutOff,numtestIts+1) ##+1 to actually get that number of testIts.
#		results, tempResult1, tempResult2 = [0], 0.0, 0.0
#		h = testValues[0] - testValues[1]
#		for i in range(numtestIts):
#			tempResult2 = tempResult1 + 0.5*h*(calc_test_qqBar_int_func(testValues[i+1],S123) + calc_test_qqBar_int_func(testValues[i],S123))
#			results.append(1.0 - math.exp((-1.0)*tempResult2))
#			tempResult1 = tempResult2
#		return results
#
#	def calc_test_qg_int_func(PperpSquared,S123):
#		"""A function to return a value for the function in the Pperp^2 test integral, for qg."""
#		assert assertions.all_are_numbers([PperpSquared,S123])
#		alphaSNow = alphaS.calculate(PperpSquared)
#		term1 = (3.0*alphaSNow)/(4.0*math.pi*PperpSquared)
#		part1 = (2.0*math.log(S123/PperpSquared)) + ((5.0*PperpSquared)/S123) - ((2.0*PperpSquared*PperpSquared)/(S123*S123))
#		part2 = (PperpSquared*PperpSquared*PperpSquared/(3.0*S123*S123*S123)) - (10.0/3.0)
#		term2 = part1 + part2
#		return term1*term2
#
#	def calc_test_qg_probs(S123,numtestIts):
#		"""A function to calculate the test probabilities as a function of Pperp^2, for qg."""
#		assert assertions.all_are_numbers([S123,numtestIts])
#		intCutOff = constants.cut_off_energy()*constants.cut_off_energy()
#		topLimit = S123 ##PperpSquaredMax = S123 for first dipole emission.
#		testValues = numpy.linspace(topLimit,intCutOff,numtestIts+1) ##+1 to actually get that number of testIts.
#		results, tempResult1, tempResult2 = [0], 0.0, 0.0
#		h = testValues[0] - testValues[1]
#		for i in range(numtestIts):
#			tempResult2 = tempResult1 + 0.5*h*(calc_test_qg_int_func(testValues[i+1],S123) + calc_test_qg_int_func(testValues[i],S123))
#			results.append(1.0 - math.exp((-1.0)*tempResult2))
#			tempResult1 = tempResult2
#		return results
#
#	def calc_test_gg_int_func(PperpSquared,S123):
#		"""A function to return a value for the function in the Pperp^2 test integral, for gg."""
#		assert assertions.all_are_numbers([PperpSquared,S123])
#		alphaSNow = alphaS.calculate(PperpSquared)
#		term1 = (3.0*alphaSNow)/(4.0*math.pi*PperpSquared)
#		part1 = (2.0*math.log(S123/PperpSquared)) + ((6.0*PperpSquared)/S123) - ((3.0*PperpSquared*PperpSquared)/(S123*S123))
#		part2 = (2.0*PperpSquared*PperpSquared*PperpSquared/(3.0*S123*S123*S123)) - (11.0/3.0)
#		term2 = part1 + part2
#		return term1*term2
#
#	def calc_test_gg_probs(S123,numtestIts):
#		"""A function to calculate the test probabilities as a function of Pperp^2, for gg."""
#		assert assertions.all_are_numbers([S123,numtestIts])
#		intCutOff = constants.cut_off_energy()*constants.cut_off_energy()
#		topLimit = S123 ##PperpSquaredMax = S123 for first dipole emission.
#		testValues = numpy.linspace(topLimit,intCutOff,numtestIts+1) ##+1 to actually get that number of testIts.
#		results, tempResult1, tempResult2 = [0], 0.0, 0.0
#		h = testValues[0] - testValues[1]
#		for i in range(numtestIts):
#			tempResult2 = tempResult1 + 0.5*h*(calc_test_gg_int_func(testValues[i+1],S123) + calc_test_gg_int_func(testValues[i],S123))
#			results.append(1.0 - math.exp((-1.0)*tempResult2))
#			tempResult1 = tempResult2
#		return results
#
#	##Begin testing:##
#	print "\n----------------------------------------------------------------------"
#	print "----------------------------------------------------------------------\n"
#	print "///////////////////////"
#	print "Testing sudakov module:"
#	print "///////////////////////"
#	assertions.pause(__name__)
#	
#	###Setup here:##
#	##Both the three test probs and allPperpsSquared run from larger Pperpsquared to smaller.
#	print "\nGenerating test values..."
#	##Generate splitting results.##
#	cutOff = constants.cut_off_energy()*constants.cut_off_energy()
#	testS123 = 63.0
#	testY = 0.756425
#	qqBarPperpSquareds, qqBarYs = [], []
#	qgPperpSquareds, qgYs = [], []
#	ggPperpSquareds, ggYs = [], []
#	testIts = 100000.0
#	numberBins = 100
#	for n in range(int(testIts)):
#		var1, var2 = solve_for_qqBar(testS123,testS123)
#		qqBarPperpSquareds.append(var1)
#		qqBarYs.append(var2)
#		var3, var4 = solve_for_qg(testS123,testS123)
#		qgPperpSquareds.append(var3)
#		qgYs.append(var4)
#		var5, var6 = solve_for_gg(testS123,testS123)
#		ggPperpSquareds.append(var5)
#		ggYs.append(var6)
#	qqBarFellBelow = qqBarPperpSquareds.count(None)
#	qqBarPperpSquareds = [x for x in qqBarPperpSquareds if x is not None]
#	qgFellBelow = qgPperpSquareds.count(None)
#	qgPperpSquareds = [x for x in qgPperpSquareds if x is not None]
#	ggFellBelow = ggPperpSquareds.count(None)
#	ggPperpSquareds = [x for x in ggPperpSquareds if x is not None]
#	bins = numpy.linspace(cutOff,testS123,numberBins+1)
#	qqBarPperpSquaredsHist, qqBarBins = numpy.histogram(qqBarPperpSquareds,bins)
#	qgPperpSquaredsHist, qgBins = numpy.histogram(qgPperpSquareds,bins)
#	ggPperpSquaredsHist, ggBins = numpy.histogram(ggPperpSquareds,bins)
#	qqBarPperpSquareds2, qgPperpSquareds2, ggPperpSquareds2 = [], [], []
#	for i in range(numberBins):
#		qqBarPperpSquareds2.insert(0,sum(qqBarPperpSquaredsHist[-(i+1):])/testIts) ##Use insert at 0 to keep the reverse ordering.
#		qgPperpSquareds2.insert(0,sum(qgPperpSquaredsHist[-(i+1):])/testIts)
#		ggPperpSquareds2.insert(0,sum(ggPperpSquaredsHist[-(i+1):])/testIts)
#	allGenerated = [qqBarPperpSquareds2,qgPperpSquareds2,ggPperpSquareds2]
#	allBins = [qqBarBins[:-1],qgBins[:-1],ggBins[:-1]] ##As allGenerated has the number created up to the end of the bin.
#	print "\nqqBar probability is", qqBarPperpSquareds2[0]
#	print "Adding the", qqBarFellBelow, "that fell below the cutoff gives a final probability of", qqBarPperpSquareds2[0] + (qqBarFellBelow/testIts)
#	if precision.check_numbers_equal(qqBarPperpSquareds2[0] + (qqBarFellBelow/testIts),1.0):
#		print "Test completed successfully!"
#	else:
#		print "Not consistent, test failed!"
#	print "\nqg probability is", qgPperpSquareds2[0]
#	print "Adding the", qgFellBelow, "that fell below the cutoff gives a final probability of", qgPperpSquareds2[0] + (qgFellBelow/testIts)
#	if precision.check_numbers_equal(qgPperpSquareds2[0] + (qgFellBelow/testIts),1.0):
#		print "Test completed successfully!"
#	else:
#		print "Not consistent, test failed!"
#	print "\ngg probability is", ggPperpSquareds2[0]
#	print "Adding the", ggFellBelow, "that fell below the cutoff gives a final probability of", ggPperpSquareds2[0] + (ggFellBelow/testIts)
#	if precision.check_numbers_equal(ggPperpSquareds2[0] + (ggFellBelow/testIts),1.0):
#		print "Test completed successfully!"
#	else:
#		print "Not consistent, test failed!"
#	assertions.pause(__name__)
#	##Generate prob functions.##
#	##Using testIts from above as the number of integration slices.
#	allPperpsSquared = numpy.linspace(testS123,cutOff,testIts+1) ##+1 to actually get that number of testIts.
#	qqBarProbs, qgProbs = calc_test_qqBar_probs(testS123,int(testIts)), calc_test_qg_probs(testS123,int(testIts))
#	ggProbs = calc_test_gg_probs(testS123,int(testIts))
#	##Don't forget the test values are the probability of a decay by a value, not at that given value.
#	qqBarProbsMax, qgProbsMax, ggProbsMax = qqBarProbs[-1], qgProbs[-1], ggProbs[-1]
#	allProbs = [qqBarProbs,qgProbs,ggProbs]
#	testParticleCodes = [[1,-1],[2,-2],[3,-3],[4,-4],[5,-5],[6,-6],[-1,1],[-2,2],[-3,3],[-4,4],[-5,5],[-6,6],[21,21]]
#	testParticleCodes += [[1,21],[2,21],[3,21],[4,21],[5,21],[6,21],[-1,21],[-2,21],[-3,21],[-4,21],[-5,21],[-6,21]]
#	testParticleCodes += [[21,1],[21,2],[21,3],[21,4],[21,5],[21,6],[21,-1],[21,-2],[21,-3],[21,-4],[21,-5],[21,-6]]
#
#	##Test completeness. i.e all split or fall below:##
#	print "\n--------------------------------------------------\n"
#	print "Testing completeness. i.e all split or fall below:\n"
#	print "\nFor qqBar,", qqBarFellBelow, "fell below the cutoff energy."
#	print len(qqBarPperpSquareds), "successfully split."
#	if ((testIts - len(qqBarPperpSquareds)) == qqBarFellBelow):
#		print "Consistent, test completed successfully!"
#	else:
#		print "Not consistent, test failed!"
#	print "\nFor qg,", qgFellBelow, "fell below the cutoff energy."
#	print len(qgPperpSquareds), "successfully split."
#	if ((testIts - len(qgPperpSquareds)) == qgFellBelow):
#		print "Consistent, test completed successfully!"
#	else:
#		print "Not consistent, test failed!"
#	print "\nFor gg,", ggFellBelow, "fell below the cutoff energy."
#	print len(ggPperpSquareds), "successfully split."
#	if ((testIts - len(ggPperpSquareds)) == ggFellBelow):
#		print "Consistent, test completed successfully!"
#	else:
#		print "Not consistent, test failed!"
#	print "\nFinished testing completeness. i.e all split or fall below."
#	assertions.pause(__name__)
#
#	##Test individual sudakovs:##
#	print "\n--------------------------------------------------\n"
#	print "Testing individual sudakovs:\n"
#	xaxes = [r"$P_{\perp}^{2}$",r"$P_{\perp}^{2}$",r"$P_{\perp}^{2}$"]
#	yaxes = [r"P(split by $P_{\perp}^{2}$)",r"P(split by $P_{\perp}^{2}$)",r"P(split by $P_{\perp}^{2}$)"]
#	titles = ['qqBar','qg','gg']
#	figure,axes = pyplot.subplots(2,2)
#	axes = axes.ravel()
#	for idx,ax in enumerate(axes):
#		if (idx == 3):
#			continue
#		#ax.hist(histData[idx], numberBins, weights = histWeights[idx], alpha = 0.5) #normed = 1
#		ax.plot(allBins[idx], allGenerated[idx], linestyle = "solid", color = "blue", linewidth = 2)
#		ax.plot(allPperpsSquared, allProbs[idx], linestyle = "solid", color = "red", linewidth = 2)
#		ax.set_title(titles[idx])
#		ax.set_xlabel(xaxes[idx])
#		ax.set_ylabel(yaxes[idx])
#		ax.axvline(cutOff, linestyle = '--', color = 'red')
#	pyplot.suptitle(r"Generated probabilities of a split by $P_{\perp}^{2}$ vs expected values for:" + str(int(testIts)) + " iterations")
#	line1 = pyplot.Line2D((0,1),(0,0), color="blue", linewidth = 2)
#	line2 = pyplot.Line2D((0,1),(0,0), color="red", linewidth = 2)
#	pyplot.figlegend([line1,line2],["Generated","Expected"],"lower right")
#	figure.delaxes(axes[3]) ##Remove blank fourth subplot.
#	pyplot.tight_layout()
#	pyplot.subplots_adjust(top=0.85) #Space main title out to prevent overlapping.
#	assertions.show_graph()
#	print "\nFinished testing individual sudakovs."
#	assertions.pause(__name__)
#
#	##Test solve_g_emission and solve:##
#	print "\n--------------------------------------------------\n"
#	print "Testing solve_g_emission and solve:\n"
#	get_name = particleData.knownParticles.get_name_from_code
#	print "Using S123 and PperpSquaredMax =", testS123
#	counter = 0
#	for testCodes in testParticleCodes:
#		counter += 1
#		print "\nCalling solve_g_emission with", get_name(testCodes[0]), "and", get_name(testCodes[1]), "gives:"
#		print solve_g_emission(testS123,testS123,testCodes[0],testCodes[1])
#		print "Calling solve with", get_name(testCodes[0]), "and", get_name(testCodes[1]), "gives:"
#		print solve(testS123,testS123,testCodes[0],testCodes[1])
#		print "--------------------"
#		if (counter%2 == 0):
#			assertions.pause(__name__)
#	print "\nFinished testing solve_g_emission and solve."
#	assertions.pause(__name__)
#
#	##Test calc_G2 and calc_inv_G2:##
#	print "\n--------------------------------------------------\n"
#	print "Testing calc_G2 and calc_inv_G2:\n"
#	print "Calling calc_G2 with:", testY
#	print "returns:", calc_G2(testY)
#	print "Calling calc_inv_G2 with :", calc_G2(testY)
#	print "returns:", calc_inv_G2(calc_G2(testY))
#	if (calc_inv_G2(calc_G2(testY)) == testY):
#		print "Test successful!"
#	else:
#		print "Test failed!"
#	print "\nFinished testing calc_G2 and calc_inv_G2."
#	assertions.pause(__name__)
#
#	##All other functions tested implicitly in the above tests.

	##Done testing:##
	print "\n---------------------------------------------\n"
	print "/////////////////////////////////"
	print "Finished checking sudakov module!"
	print "/////////////////////////////////"