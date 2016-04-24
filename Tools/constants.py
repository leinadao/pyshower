####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module containing constants to be used throughout dipole showering."""

##Import required modules:##

print "\n/////////////////////////"
print "Loading constants module:"
print "/////////////////////////\n"

##Functions:##

###~~~Defining program values~~~###

def number_decimal_places():
	"""A function to define the number of decimal places being used for tests throughout."""
	__numDecPlaces = 7 ##Reduced from 11 to stop LHEF values not giving massless particles.
	return __numDecPlaces

def division_number_decimal_places():
	"""A function to define the number of decimal places being used throughout when checking for division by zero."""
	##Python precision reached: Lorentz test of Pi^2 conservation were failing so reduced to 11.
	__numDecPlaces = 11
	return __numDecPlaces

def LHEF_ME_number_decimal_places():
	"""A function to define the number of decimal places being used when saving matrix elements in LHEF format."""
	__numDecPlaces = 10 ##Using the value that currently seems to work best for Pythia.
	return __numDecPlaces

def LHEF_DS_number_decimal_places():
	"""A function to define the number of decimal places being used when saving showered results in LHEF format."""
	__numDecPlaces = 10 ##Using the value that currently seems to work best for Pythia.
	return __numDecPlaces

def Nc():
	"""A function to return the value of Nc, the number of colours, used throughout."""
	__Nc = 3.0
	return __Nc

def Tr():
	"""Return a value for the colour factor Tr used throughout."""
	##For g -> q+qBar.
	__Tr = 0.5
	return __Tr

def one_loop_alphaS_of_Mz_squared():
	"""Return a value for alphaS(Mz^2) from literature/experiment for the one loop strong coupling constant."""
	##Source: Professor Frank Krauss' Book (currently unpublished).
	__oneLoopAlphaSOfMzSquared = 0.118
	return __oneLoopAlphaSOfMzSquared

def one_loop_alphaEM_of_Mz_squared():
	"""Return a value for alphaEM(Mz^2) from literature/experiment for the one loop fine structure constant."""
	##Source: Paper 'Measurement of the Running of the Fine-Structure Constant'.
	__oneLoopAlphaEMOfMzSquared = 1.0/128.886
	return __oneLoopAlphaEMOfMzSquared

def two_loop_alphaS_of_Mz():
	"""Return a value for alphaS(Mz) from literature/experiment for the two loop strong coupling constant."""
	##Source: Paper 'Measurements of the strong coupling constant and the QCD colour
	##Paper 'Four-jet observables from hadronic Z decays', p2-5. (for factors)
	##Normally you use Mz^2 as your reference value and a function of Q^2 but here they have used Mz and Q.
	__twoLoopAlphaSOfMz = 0.118
	return __twoLoopAlphaSOfMz

def four_loop_lambda():
	"""Return a value for lambda from literature/experiment for the four loop strong coupling constant."""
	##Source: Professor Frank Krauss' Book (currently unpublished).
	__fourLoopLambda = 0.208364759205
	return __fourLoopLambda

def cut_off_energy():
	"""Return a value for the cut off energy at which showering ends and hadronisation begins."""
	##Chosen after viewing multiple sources. To be varied to fit hadronisation results.
	__cutOff = 1.0
	return __cutOff

def aperys_constant():
	"""Return a value for aperys constant."""
	##Source: Wikipedia page on 'Apery's Constant'.
	__aperysConstant = 1.20205690315959428539973816151144999076498629
	return __aperysConstant

def sin_squared_thetaW():
	"""A function to return a value for sin^2 of the weak mixing angle."""
	##Source: (5) PDG - Physical Constants PDF.
	__sSThetaW = 0.23126
	return __sSThetaW

def active_q_codes():
	"""A function to return the default active q codes."""
	##Top too heavy for Z-boson, b and c massless causes problems in Pythia.
	__activeQCodes = [1,2,3]
	return __activeQCodes

###~~~Functions of defining program values~~~###

def Cf():
	"""Return a value for the colour factor Cf using the formule Cf = (Nc^2 - 1)/(2*Nc)."""
	##For q -> q+g.
	__Nc = Nc()
	__denominator1 = (2.0*__Nc)
	return (__Nc*__Nc - 1.0)/__denominator1

def Ca():
	"""Returns a value of the colour factor Ca using the formula Ca = 2*Tr*Nc."""
	##Source: 'Three-loop beta-functions for top-Yukawa and the Higgs self-interaction in the Standard Model', p6.
	##For g -> g+g.
	__Ca = 2.0 * Tr() * Nc()
	return __Ca

def neutral_weak_couplings(lowerIndex,upperIndex):
	"""A function to calculate the neutral weak couplings for e+e- -> ffBar where f is a fermion."""
	##Source: Textbook 'Dynamics of the Standard Model', p59.
	assert ((type(lowerIndex) == str) and (lowerIndex in ["v","a"]))
	assert (type(upperIndex) == int) ##Requires PDG ID's
	__gC = particleData.knownParticles.get_code_from_name
	__sSThetaW = sin_squared_thetaW()
	__set1 = [__gC("electron"),__gC("muon"),__gC("tau")]
	__set2 = [__gC("u-quark"),__gC("c-quark"),__gC("t-quark")]
	__set3 = [__gC("d-quark"),__gC("s-quark"),__gC("b-quark")]
	__set4 = [__gC("electron-neutrino"),__gC("muon-neutrino"),__gC("tau-neutrino")]
	if (lowerIndex == "a"):
		if (upperIndex in __set1):
			__nWCResult = -0.5
		elif (upperIndex in __set2):
			__nWCResult = 0.5
		elif (upperIndex in __set3):
			__nWCResult = -0.5
		elif (upperIndex in __set4):
			__nWCResult = 0.5
	elif (lowerIndex == "v"):
		if (upperIndex in __set1):
			__nWCResult = (2.0*__sSThetaW) - 0.5
		elif (upperIndex in __set2):
			__nWCResult = 0.5 - ((4.0/3.0)*__sSThetaW)
		elif (upperIndex in __set3):
			__nWCResult = ((2.0/3.0)*__sSThetaW) - 0.5
		elif (upperIndex in __set4):
			__nWCResult = 0.5
	return __nWCResult

##Classes:##

##Import here to prevent cyclic import error:
import particleData

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##
	import assertions
	import precision

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "/////////////////////////"
	print "Testing constants module:"
	print "/////////////////////////"
	assertions.pause(__name__)

	##Setup here:##
	tExpectedAlphaEM = 1.0/128.886
	tLowerIndices = ['v','a']

	##Test number_decimal_places:##
	print "\n--------------------------------------------------\n"
	print "Testing number_decimal_places:\n"
	print "Calling number_decimal_places returns: " , number_decimal_places()
	if number_decimal_places() == 7:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print " \nFinished testing number_decimal_places."
	assertions.pause(__name__)

	##Test division_number_decimal_places:##
	print "\n--------------------------------------------------\n"
	print "Testing division_number_decimal_places:\n"
	print "Calling division_number_decimal_places returns: " , division_number_decimal_places()
	if division_number_decimal_places() == 11:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print " \nFinished testing division_number_decimal_places."
	assertions.pause(__name__)

	##Test LHEF_ME_number_decimal_places:##
	print "\n--------------------------------------------------\n"
	print "Testing LHEF_ME_number_decimal_places:\n"
	print "Calling LHEF_ME_number_decimal_places returns: " , LHEF_ME_number_decimal_places()
	if LHEF_ME_number_decimal_places() == 10:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print " \nFinished testing LHEF_ME_number_decimal_places."
	assertions.pause(__name__)

	##Test LHEF_DS_number_decimal_places:##
	print "\n--------------------------------------------------\n"
	print "Testing LHEF_DS_number_decimal_places:\n"
	print "Calling LHEF_DS_number_decimal_places returns: " , LHEF_DS_number_decimal_places()
	if LHEF_DS_number_decimal_places() == 10:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print " \nFinished testing LHEF_DS_number_decimal_places."
	assertions.pause(__name__)

	##Test Nc:##
	print "\n--------------------------------------------------\n"
	print "Testing Nc:\n"
	print "Nc = " , Nc()
	if Nc() == 3:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing Nc."
	assertions.pause(__name__)

	##Test Tr:##
	print "\n--------------------------------------------------\n"
	print "Testing Tr:\n"
	print "Tr = " , Tr()
	if Tr() == 0.5:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing Tr."
	assertions.pause(__name__)

	##Test one_loop_alphaS_of_Mz_squared:##
	print "\n--------------------------------------------------\n"
	print "Testing one_loop_alphaS_of_Mz_squared:\n"
	print "one_loop_alphaS_of_Mz_squared = " , one_loop_alphaS_of_Mz_squared()
	if one_loop_alphaS_of_Mz_squared() == 0.118:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing one_loop_alphaS_of_Mz_squared."
	assertions.pause(__name__)

	##Test one_loop_alphaEM_of_Mz_squared:##
	print "\n--------------------------------------------------\n"
	print "Testing one_loop_alphaEM_of_Mz_squared:\n"
	print "one_loop_alphaEM_of_Mz_squared = " , one_loop_alphaEM_of_Mz_squared()
	if one_loop_alphaEM_of_Mz_squared() == (tExpectedAlphaEM):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing one_loop_alphaEM_of_Mz_squared."
	assertions.pause(__name__)

	##Test two_loop_alphaS_of_Mz:##
	print "\n--------------------------------------------------\n"
	print "Testing two_loop_alphaS_of_Mz:\n"
	print "two_loop_alphaS_of_Mz = " , two_loop_alphaS_of_Mz()
	if two_loop_alphaS_of_Mz() == 0.118:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing two_loop_alphaS_of_Mz."
	assertions.pause(__name__)

	##Test four_loop_lambda:##
	print "\n--------------------------------------------------\n"
	print "Testing four_loop_lambda:\n"
	print "four_loop_lambda = " , four_loop_lambda()
	if four_loop_lambda() == 0.208364759205:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing four_loop_lambda."
	assertions.pause(__name__)

	##Test cut_off_energy:##
	print "\n--------------------------------------------------\n"
	print "Testing cut_off_energy:\n"
	print "cut_off_energy() returns:" , cut_off_energy()
	if cut_off_energy() == 1.0:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing cut_off_energy."
	assertions.pause(__name__)

	##Test aperys_constant:##
	print "\n--------------------------------------------------\n"
	print "Testing aperys_constant:\n"
	print "aperys_constant() returns:" , aperys_constant()
	if aperys_constant() == 1.20205690315959428539973816151144999076498629:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing aperys_constant."
	assertions.pause(__name__)

	##Test sin_squared_thetaW:##
	print "\n--------------------------------------------------\n"
	print "Testing sin_squared_thetaW:\n"
	print "sin_squared_thetaW() returns:" , sin_squared_thetaW()
	if sin_squared_thetaW() == 0.23126:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing sin_squared_thetaW."
	assertions.pause(__name__)

	##Test active_q_codes:##
	print "\n--------------------------------------------------\n"
	print "Testing active_q_codes:\n"
	print "active_q_codes() returns:" , active_q_codes()
	if active_q_codes() == [1,2,3]:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing active_q_codes."
	assertions.pause(__name__)

	##Test Cf:##
	print "\n--------------------------------------------------\n"
	print "Testing Cf:\n"
	print "Nc = " , Nc()
	print "Cf = " , Cf()
	if Cf() == ((Nc()*Nc()-1.0)/(2.0*Nc())):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing Cf."
	assertions.pause(__name__)

	##Test Ca:##
	print "\n--------------------------------------------------\n"
	print "Testing Ca:\n"
	print "Nc = " , Nc()
	print "Ca = " , Ca()
	if Ca() == 2.0*Tr()*Nc():
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing Ca."
	assertions.pause(__name__)

	##Test neutral_weak_couplings:##
	print "\n--------------------------------------------------\n"
	print "Testing neutral_weak_couplings:\n"
	tFermionCodes = particleData.knownParticles.get_known_fermions()
	tG_N = particleData.knownParticles.get_name_from_code
	results = []
	for fC in tFermionCodes:
		for lI in tLowerIndices:
			print "Calling with lower index", lI ,"and upper index", tG_N(fC), "returns:"
			print neutral_weak_couplings(lI,fC)
			results.append(neutral_weak_couplings(lI,fC))
			assertions.pause(__name__)
	tExpected = [-0.34582666666,-0.5,0.19165333333,0.5,-0.34582666666,-0.5,0.19165333333,0.5,-0.34582666666,-0.5,0.19165333333,0.5]
	tExpected += [-0.03748,-0.5,0.5,0.5,-0.03748,-0.5,0.5,0.5,-0.03748,-0.5,0.5,0.5]
	tSuccessful = True
	for i, v in enumerate(results):
		if not precision.check_numbers_equal(results[i],tExpected[i]):
			tSuccessful = False
	if tSuccessful:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing neutral_weak_couplings."
	assertions.pause(__name__)
	
	##Done testing:##
	print "\n---------------------------------------------\n"
	print "///////////////////////////////////"
	print "Finished checking constants module!"
	print "///////////////////////////////////"