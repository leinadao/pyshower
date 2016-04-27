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
import dataLoggers
import particleData
import runningCouplings
import kinematics

print "\n///////////////////////"
print "Loading sudakov module:"
print "///////////////////////\n"

##Initiate coupling constants:##
alphaS = runningCouplings.oneLoopAlphaS()
alphaEM = runningCouplings.oneLoopAlphaEM()

##Initialise data loggers:##
crossSecsGluSplit = dataLoggers.dataLogger()
crossSecsGluProd = dataLoggers.dataLogger()

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

def solve_sudakovs(PperpSquaredMax,S123,difCrossSecIds,code1,code2,logCrosSecs = None):
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
	if (logCrosSecs == 1):
		crossSecsGluSplit.store(qg_to_qQQBar_cS(S123,__PperpSquaredi,__yi))
		crossSecsGluProd.store(qg_to_qgg_cS(S123,__PperpSquaredi,__yi))
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
	__PperpSquared, __y, __difCrossSecId = solve_sudakovs(PperpSquaredMax,S123,__difCrossSecIds,code1,code2,1)
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

	##Module re-written so old tests removed.##
	print "Module re-written so old tests removed."

	##Done testing:##
	print "\n---------------------------------------------\n"
	print "/////////////////////////////////"
	print "Finished checking sudakov module!"
	print "/////////////////////////////////"