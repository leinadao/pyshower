####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module for handling Lorentz boosts and rotations of four-vectors."""

##Space-like vectors (i.e P.P = M^2 < 0) are not treated as not needed for a final state shower.

##Import required modules:##
import math
import assertions
import precision
import fourVectors

print "\n///////////////////////"
print "Loading lorentz module:"
print "///////////////////////\n"

##Functions:##

def check_is_lorentzBoost(toCheck):
	"""A function to check for an instance of the lorentzBoost class."""
	return isinstance(toCheck,lorentzBoost)

def check_is_lorentzRotation(toCheck):
	"""A function to check for an instance of the lorentzRotation class."""
	return isinstance(toCheck,lorentzRotation)

def check_is_boostAndRotate(toCheck):
	"""A function to check for an instance of the boostAndRotate class."""
	return isinstance(toCheck,boostAndRotate)

def check_timelike(pV):
	"""A function to check whether a momentum vector is time-like (i.e P.P = M^2 > 0)."""
	assert fourVectors.check_is_fourVector(pV)
	assert pV.__nonzero__()
	__MinkowskiSP = pV * pV
	__MinkowskiSP = precision.round_python_errors(__MinkowskiSP)
	return (__MinkowskiSP > 0.0)

def check_lightlike(pV):
	"""A function to check whether a momentum vector is light-like (i.e P.P = M^2 = 0)."""
	assert fourVectors.check_is_fourVector(pV)
	assert pV.__nonzero__()
	__MinkowskiSP = pV * pV
	return precision.check_numbers_equal(__MinkowskiSP,0.0)

def calculate_m_from_p_vector(pV):
	"""A function to calculate the mass of a particle from its momentum four-vector."""
	assert fourVectors.check_is_fourVector(pV)
	assert pV.__nonzero__()
	return math.sqrt(pV * pV)
	
def check_v_not_faster_than_c(scalarVelocity):
	"""A function to check that a scalar velocity isn't faster than the speed of light."""
	##Allows equal to the speed of light to accommodate massless particles.
	assert assertions.all_are_numbers([scalarVelocity])
	__scalarVelocity = assertions.force_float_number(scalarVelocity)
	return (__scalarVelocity <= 1.0)

def p_vector_to_v(pV):
	"""A function to return the velocity four-vector from a momentum four-vector."""
	assert fourVectors.check_is_fourVector(pV)
	assert pV.__nonzero__()
	assert check_timelike(pV)
	__mass = calculate_m_from_p_vector(pV)
	__velocityResult = pV.copy()
	__velocityResult /= __mass
	return __velocityResult

##Classes:##

class lorentzBoost(object):
	"""A class for performing Lorentz boosts on fourVector objects."""

	def __init__(self,momentumOfCMF):
		"""A function to initiate a Lorentz boost to and from a centre of mass frame (E,p)."""
		assert fourVectors.check_is_fourVector(momentumOfCMF)
		assert momentumOfCMF.__nonzero__()
		assert (check_timelike(momentumOfCMF) or check_lightlike(momentumOfCMF))
		self.__momentumOfCMF = momentumOfCMF.copy()
		self.__reverseMomentumOfCMF = momentumOfCMF.copy()
		for __i1 in range(3):
			self.__reverseMomentumOfCMF[__i1+1] *= (-1.0)
		self.__massOfCMF = calculate_m_from_p_vector(momentumOfCMF)

	def get_CMF_mass(self):
		"""A function to return the mass associated with the CMF, in the lab frame."""
		return self.__massOfCMF

	def get_beta(self):
		"""A function to return the beta vector associated with the centre of mass frame."""
		if check_timelike(self.__momentumOfCMF):
			__beta = p_vector_to_v(self.__momentumOfCMF)
		elif check_lightlike(self.__momentumOfCMF):
			print "Momentum Vector is light-like. Can't handle."
			assert check_timelike(self.__momentumOfCMF)
		else:
			print "Momentum Vector is space-like. Can't handle."
			assert check_timelike(self.__momentumOfCMF)
		return __beta

	def get_momentum_in_CMF(self):
		"""A function to return the four-vector in its COM frame."""
		return fourVectors.fourVector(self.__massOfCMF, 0.0, 0.0, 0.0)

	def boost(self,pVIn):
		"""A function to perform the lorentz boost on a given (E,p) four-vector."""
		assert fourVectors.check_is_fourVector(pVIn)
		assert pVIn.__nonzero__()
		##Make sure it's not space-like:
		__pVIn = pVIn.copy()
		assert (check_timelike(__pVIn) or check_lightlike(__pVIn))
		##Make sure not trying to boost into its own CMF:
		if check_lightlike(__pVIn):
			assert (__pVIn != self.__momentumOfCMF)
		__pVOut = fourVectors.fourVector()
		__energyPrime = (self.__momentumOfCMF * __pVIn)/self.__massOfCMF
		__alpha = (__pVIn[0]+__energyPrime)/(self.__massOfCMF + self.__momentumOfCMF[0])
		__pVOut[0] = __energyPrime
		for __i2 in range(3):
			__pVOut[__i2+1] = __pVIn[__i2+1] - __alpha * self.__momentumOfCMF[__i2+1]
		##Make sure its still not space-like after the boost:
		assert (check_timelike(__pVOut) or check_lightlike(__pVOut))
		return __pVOut

	def inverse_boost(self,pVIn):
		"""A function to perform the inverse lorentz boost on a given (E,p) four-vector."""
		assert fourVectors.check_is_fourVector(pVIn)
		assert pVIn.__nonzero__()
		##Make sure it's not space-like:
		assert (check_timelike(pVIn) or check_lightlike(pVIn))
		__pVOut = fourVectors.fourVector()
		__energyOut = (self.__reverseMomentumOfCMF * pVIn)/self.__massOfCMF
		__alpha = (pVIn[0]+__energyOut)/(self.__massOfCMF + self.__momentumOfCMF[0])
		__pVOut[0] = __energyOut
		for __i3 in range(3):
			__pVOut[__i3+1] = pVIn[__i3+1] - __alpha * self.__reverseMomentumOfCMF[__i3+1]
		##Make sure its still not space-like after the boost:
		assert (check_timelike(__pVOut) or check_lightlike(__pVOut))
		return __pVOut

	def __mul__(self,pVToBoost):
		"""A function to call a lorentz boost on the four-vector multiplied by."""
		assert fourVectors.check_is_fourVector(pVToBoost)
		assert pVToBoost.__nonzero__()
		return self.boost(pVToBoost)

	def __div__(self,pVToInverseBoost):
		"""A function to call an inverse lorentz boost on the four-vector divided by."""
		assert fourVectors.check_is_fourVector(pVToInverseBoost)
		assert pVToInverseBoost.__nonzero__()
		return self.inverse_boost(pVToInverseBoost)

class lorentzRotation(object):
	"""A class for performing lorentz rotations on four-vectors."""

	def calculate_angles(self):
		"""A function to calculate the angles required to rotate the initialising vector onto the z-axis"""
		##The factor of -1 is to rotate onto the z-axis given that the angles are those from it / the x-axis:
		self.__theta = -1.0 * math.acos(self.__initialisingVector[3]/self.__r)
		##atan2 is calculating: atan(self.__initialisingVector[2]/self.__initialisingVector[1]):
		self.__phi = math.atan2(self.__initialisingVector[2],self.__initialisingVector[1])
		if self.__phi < 0.0:
			self.__phi += 2.0*math.pi
		self.__phi *= -1.0

	def __init__(self, initialisingVector):
		"""A function to initialise the lorentzRotation class."""
		assert fourVectors.check_is_fourVector(initialisingVector)
		assert initialisingVector.__nonzero__()
		self.__initialisingVector = initialisingVector.copy()
		self.__r = self.__initialisingVector.calculate_cartesian_magnitude()
		self.calculate_angles()

	def rotate(self,vectorToRotate):
		"""A function to apply the forward rotation to the given four-vector."""
		##This is the rotation initialising four-vector onto z-axis.
		assert fourVectors.check_is_fourVector(vectorToRotate)
		assert vectorToRotate.__nonzero__()
		__v2R, __c, __s, __t, __p, = vectorToRotate.copy(), math.cos , math.sin, self.__theta, self.__phi
		__rotatedResult = vectorToRotate.copy()
		__rotatedResult[1] = __c(__t)*__c(__p)*__v2R[1] - __c(__t)*__s(__p)*__v2R[2] + __s(__t)*__v2R[3]
		__rotatedResult[2] = __s(__p)*__v2R[1] + __c(__p)*__v2R[2]
		__rotatedResult[3] = -__s(__t)*__c(__p)*__v2R[1] + __s(__t)*__s(__p)*__v2R[2] + __c(__t)*__v2R[3]
		return __rotatedResult

	def inverse_rotate(self,vectorToInverseRotate):
		"""A function to apply the backwards rotation to the given four-vector."""
		##This is the rotation initialising four-vector from the z-axis.
		assert fourVectors.check_is_fourVector(vectorToInverseRotate)
		assert vectorToInverseRotate.__nonzero__()
		__v2R, __c, __s, __t, __p, = vectorToInverseRotate.copy(), math.cos , math.sin, self.__theta, self.__phi
		__rotatedResult = vectorToInverseRotate.copy()
		__rotatedResult[1] = __c(__p)*__c(__t)*__v2R[1] + __s(__p)*__v2R[2] - __c(__p)*__s(__t)*__v2R[3]
		__rotatedResult[2] = -__s(__p)*__c(__t)*__v2R[1] + __c(__p)*__v2R[2] + __s(__p)*__s(__t)*__v2R[3]
		__rotatedResult[3] = __s(__t)*__v2R[1] + __c(__t)*__v2R[3]
		return __rotatedResult

	def __mul__(self,vectorToRotate):
		"""A function to call a lorentz rotation on the four-vector multiplied by."""
		assert fourVectors.check_is_fourVector(vectorToRotate)
		assert vectorToRotate.__nonzero__()
		return self.rotate(vectorToRotate)

	def __div__(self,vectorToInverseRotate):
		"""A function to call an inverse lorentz rotation on the four-vector divided by."""
		assert fourVectors.check_is_fourVector(vectorToInverseRotate)
		assert vectorToInverseRotate.__nonzero__()
		return self.inverse_rotate(vectorToInverseRotate)

class boostAndRotate(object):
	"""A function to perform a combined boost and rotation for two particles to and from their orientated CMF."""

	def __init__(self,pV1,pV2,numToRecoil):
		"""A function to initialise the joint boost and rotation."""
		assert fourVectors.check_is_fourVector(pV1)
		assert fourVectors.check_is_fourVector(pV2)
		assert (pV1.__nonzero__() and pV2.__nonzero__())
		assert assertions.valid_dipole_index(numToRecoil) ##Use as the same convention as dipoles of 0,1.
		self.__vector1, self.__vector2 = pV1.copy(), pV2.copy()
		self.__momentumOfCMF = self.__vector1 + self.__vector2
		self.__theBoost = lorentzBoost(self.__momentumOfCMF)
		##numToRecoil used to decide which to rotate onto the -ve z-axis (i.e recoiling).
		if (numToRecoil == 0): ##Opposite of numToRecoil used, as to be on +ve z-axis
			self.__boostedVector = self.__theBoost*self.__vector2
			self.__theRotation = lorentzRotation(self.__boostedVector)
		else:
			self.__boostedVector = self.__theBoost*self.__vector1
			self.__theRotation = lorentzRotation(self.__boostedVector)

	def boost_rotate(self,vectorIn):
		"""A function to apply the forward boost and rotation to the given four-vector."""
		assert fourVectors.check_is_fourVector(vectorIn)
		assert vectorIn.__nonzero__()
		__vectorIn = vectorIn.copy()
		return self.__theRotation*(self.__theBoost*__vectorIn)

	def inverse_boost_rotate(self,vectorOut):
		"""A function to apply the backwards boost and rotation to the given four-vector."""
		assert fourVectors.check_is_fourVector(vectorOut)
		assert vectorOut.__nonzero__()
		__vectorOut = vectorOut.copy()
		return self.__theBoost / (self.__theRotation / __vectorOut)

	def __mul__(self,vectorIn):
		"""A function to call a lorentz boost and rotation on the four-vector multiplied by."""
		assert fourVectors.check_is_fourVector(vectorIn)
		assert vectorIn.__nonzero__()
		__vectorIn = vectorIn.copy()
		return self.boost_rotate(__vectorIn)

	def __div__(self,vectorOut):
		"""A function to call an inverse lorentz boost and rotation on the four-vector divided by."""
		assert fourVectors.check_is_fourVector(vectorOut)
		assert vectorOut.__nonzero__()
		__vectorOut = vectorOut.copy()
		return self.inverse_boost_rotate(__vectorOut)

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##
	import random

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "///////////////////////"
	print "Testing lorentz module:"
	print "///////////////////////"
	assertions.pause(__name__)
	
	##Setup here:##
	print "\nGenerating test values..."
	tLightlike = fourVectors.fourVector(math.sqrt(3*3 + 4*4 + 5*5),3,4,5)
	tMomentum4V1 = fourVectors.fourVector(200,30,50,40)
	tMomentum4V2 = fourVectors.fourVector(215,66,12,52)
	tNotTooFast1 = 0.87
	tNotTooFast2 = 1.0
	tTooFast1 = 1.01
	tTooFast2 = 109
	tBoost1 = lorentzBoost(tMomentum4V1)
	tFourVectors = {'a':None,'b':None,'c':None,'d':None}
	tMomentumSquares = {'a':None,'b':None,'c':None,'d':None}
	tBoostedMomentumSquares = {'a':None,'b':None,'c':None,'d':None}
	tSumMomentums = fourVectors.fourVector(0,0,0,0)
	tBoostedSumMomentums = fourVectors.fourVector(0,0,0,0)
	##Generate random test four-vectors:##
	for tFourVector in tFourVectors:
		tX0 = random.randrange(1000,2000)/10.0
		tX1 = random.randrange(0,200)/10.0
		tX2 = random.randrange(0,200)/10.0
		tX3 = random.randrange(0,200)/10.0
		tFourVectors[tFourVector] = fourVectors.fourVector(tX0,tX1,tX2,tX3)
		tSumMomentums += tFourVectors[tFourVector]
	tBoost3 = lorentzBoost(tSumMomentums)
	for tFourVector in tFourVectors:
		tMomentumSquares[tFourVector] = tFourVectors[tFourVector] * tFourVectors[tFourVector]
		tBoostedMomentumSquares[tFourVector] = (tBoost3 * tFourVectors[tFourVector]) * (tBoost3 * tFourVectors[tFourVector])
		tSumMomentums += tFourVectors[tFourVector]
		tBoostedSumMomentums +=  (tBoost3 * tFourVectors[tFourVector])
	tRotate1 = lorentzRotation(tMomentum4V1)
	tMomentum4V1Reversed = tMomentum4V1.copy()
	for ti in range(3):
		tMomentum4V1Reversed[ti + 1] *= -1.0
	tTogetherCOM = tFourVectors['a'].copy() + tFourVectors['b']
	tBoost4 = lorentzBoost(tTogetherCOM)
	tRotate2 = lorentzRotation(tBoost4*tFourVectors['a'])
	tTogether2V1 = tMomentum4V1.copy()
	tTogether2V2 =  fourVectors.fourVector(math.sqrt(tTogether2V1*tTogether2V1 + 20*20 + 40*40 + 30*30),20,40,30)
	tTogether2V2C = tTogether2V2.copy()
	tTogether2COM = tTogether2V1 + tTogether2V2
	tBoost5 = lorentzBoost(tTogether2COM)
	tRotate3 = lorentzRotation(tBoost5*tTogether2V1)
	tMasslessEnergy = math.sqrt(0.2**2 + 0.3**2 + 0.4**2)
	tMassless4V1 = fourVectors.fourVector(tMasslessEnergy,0.2,0.3,0.4)
	tMassless4V2 = fourVectors.fourVector(tMasslessEnergy,0.3,0.4,0.2)
	tMasslessCOM = tMassless4V1 + tMassless4V2
	tMasslessBoost1 = lorentzBoost(tMasslessCOM)
	tTogetherRV1Massless = tTogether2V1.copy()
	tTogetherRV2Massless = tTogether2V2.copy()
	tTogetherRV1Massless[0] = tTogetherRV1Massless.calculate_cartesian_magnitude()
	tTogetherRV2Massless[0] = tTogetherRV2Massless.calculate_cartesian_magnitude()
	tTogetherCOMMassless = tTogetherRV1Massless + tTogetherRV2Massless
	tBoost6 = lorentzBoost(tTogetherCOMMassless)
	tRotate4 = lorentzRotation(tBoost6*tTogetherRV1Massless)
	tBoostAndRotate0 = boostAndRotate(tTogetherRV1Massless,tTogetherRV2Massless,0)
	tBoostAndRotate1 = boostAndRotate(tTogetherRV1Massless,tTogetherRV2Massless,1)
	tSpacelike = fourVectors.fourVector(3.0,1.0,3.0,2.0)
	c_n_e = precision.check_numbers_equal

	##Test check_is_lorentz'':##
	print "\n--------------------------------------------------\n"
	print "Testing check_is_lorentz'':\n"
	tSuccessful = True
	print "Calling check_is_lorentzBoost on instance:" , check_is_lorentzBoost(tBoost1)
	if not check_is_lorentzBoost(tBoost1):
		tSuccessful = False
	print "Calling check_is_lorentzBoost on wrong instance:" , check_is_lorentzBoost(tRotate1)
	if check_is_lorentzBoost(tRotate1):
		tSuccessful = False
	print "Calling check_is_lorentzBoost on 1.055:" , check_is_lorentzBoost(1.055)
	if check_is_lorentzBoost(1.055):
		tSuccessful = False
	print "Calling check_is_lorentzBoost on 'word':" , check_is_lorentzBoost('word')
	if check_is_lorentzBoost('word'):
		tSuccessful = False
	print "Calling check_is_lorentzRotation on instance:" , check_is_lorentzRotation(tRotate1)
	if not check_is_lorentzRotation(tRotate1):
		tSuccessful = False
	print "Calling check_is_lorentzRotation on wrong instance:" , check_is_lorentzRotation(tBoost1)
	if check_is_lorentzRotation(tBoost1):
		tSuccessful = False
	print "Calling check_is_lorentzRotation on 1.055:" , check_is_lorentzRotation(1.055)
	if check_is_lorentzRotation(1.055):
		tSuccessful = False
	print "Calling check_is_lorentzRotation on 'word':" , check_is_lorentzRotation('word')
	if check_is_lorentzRotation('word'):
		tSuccessful = False
	print "Calling check_is_boostAndRotate on instance:" , check_is_boostAndRotate(tBoostAndRotate1)
	if not check_is_boostAndRotate(tBoostAndRotate1):
		tSuccessful = False
	print "Calling check_is_boostAndRotate on wrong instance:" , check_is_boostAndRotate(tBoost1)
	if check_is_boostAndRotate(tBoost1):
		tSuccessful = False
	print "Calling check_is_boostAndRotate on 1.055:" , check_is_boostAndRotate(1.055)
	if check_is_boostAndRotate(1.055):
		tSuccessful = False
	print "Calling check_is_boostAndRotate on 'word':" , check_is_boostAndRotate('word')
	if check_is_boostAndRotate('word'):
		tSuccessful = False
	if tSuccessful:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_is_lorentz''."
	assertions.pause(__name__)

	##Test check_timelike and check_lightlike:##
	print "\n--------------------------------------------------\n"
	print "Testing check_timelike and check_lightlike:\n"
	tSuccessful = True
	print "Calling check_lightlike on tLightlike:" , check_lightlike(tLightlike)
	if not check_lightlike(tLightlike):
		tSuccessful = False
	print "Calling check_lightlike on timelike vector:" , check_lightlike(tMomentum4V1)
	if check_lightlike(tMomentum4V1):
		tSuccessful = False
	print "Calling check_timelike on timelike vector:" , check_timelike(tMomentum4V1)
	if not check_timelike(tMomentum4V1):
		tSuccessful = False
	print "Calling check_timelike on tLightlike:" , check_timelike(tLightlike)
	if check_timelike(tLightlike):
		tSuccessful = False
	if tSuccessful:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_timelike and check_lightlike."
	assertions.pause(__name__)

	##Test calculate_m_from_p_vector, check_v_not_faster_than_c and p_vector_to_v:##
	print "\n--------------------------------------------------\n"
	print "Testing calculate_m_from_p_vector, check_v_not_faster_than_c and p_vector_to_v:\n"
	tSuccessful = True
	print "Calling calculate_m_from_p_vector on tLightlike:" , calculate_m_from_p_vector(tLightlike)
	if not precision.check_numbers_equal(calculate_m_from_p_vector(tLightlike),0.0):
		tSuccessful = False
	print "Calling check_v_not_faster_than_c on" ,  tNotTooFast1 , ":" , check_v_not_faster_than_c(tNotTooFast1)
	if not check_v_not_faster_than_c(tNotTooFast1):
		tSuccessful = False
	print "Calling check_v_not_faster_than_c on" ,  tNotTooFast2 , ":" , check_v_not_faster_than_c(tNotTooFast2)
	if not check_v_not_faster_than_c(tNotTooFast2):
		tSuccessful = False
	print "Calling check_v_not_faster_than_c on" ,  tTooFast1 , ":" , check_v_not_faster_than_c(tTooFast1)
	if check_v_not_faster_than_c(tTooFast1):
		tSuccessful = False
	print "Calling check_v_not_faster_than_c on" ,  tTooFast2 , ":" , check_v_not_faster_than_c(tTooFast2)
	if check_v_not_faster_than_c(tTooFast2):
		tSuccessful = False
	print "Using" , tMomentum4V1 , "which has mass:" , calculate_m_from_p_vector(tMomentum4V1)
	tExpected = math.sqrt((tMomentum4V1[1]*tMomentum4V1[1]) + (tMomentum4V1[2]*tMomentum4V1[2]) + (tMomentum4V1[3]*tMomentum4V1[3]))
	print "Calling p_vector_to_v gives:" , p_vector_to_v(tMomentum4V1)
	if not (p_vector_to_v(tMomentum4V1)[2] == tMomentum4V1[2] / calculate_m_from_p_vector(tMomentum4V1)):
		print "squark2"
		tSuccessful = False
	if tSuccessful:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing calculate_m_from_p_vector, check_v_not_faster_than_c and p_vector_to_v."
	assertions.pause(__name__)

	##Test lorentzBoost class:##
	print "\n--------------------------------------------------\n"
	print "///////////////////////////"
	print "Testing lorentzBoost class:"
	print "///////////////////////////"
	assertions.pause(__name__)

	##__init__() tested implicitly through the other functions.

	##Test get_CMF_mass:##
	print "\n--------------------------------------------------\n"
	print "Testing get_CMF_mass:\n"
	print "CMF momentum vector is:" , tMomentum4V1
	print "get_CMF_mass returns:" , tBoost1.get_CMF_mass()
	if (tBoost1.get_CMF_mass() == calculate_m_from_p_vector(tMomentum4V1)):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing get_CMF_mass."
	assertions.pause(__name__)

	##Test get_beta:##
	print "\n--------------------------------------------------\n"
	print "Testing get_beta:\n"
	print "CMF momentum vector is:" , tMomentum4V1
	print "get_beta returns:" , tBoost1.get_beta()
	if (tBoost1.get_beta()[2] == tMomentum4V1[2]/calculate_m_from_p_vector(tMomentum4V1)):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing get_beta."
	assertions.pause(__name__)

	##Test get_momentum_in_CMF:##
	print "\n--------------------------------------------------\n"
	print "Testing get_momentum_in_CMF:\n"
	print "CMF momentum vector is:" , tMomentum4V1
	print "get_momentum_in_CMF returns:" , tBoost1.get_momentum_in_CMF()
	tExpected = fourVectors.fourVector(calculate_m_from_p_vector(tMomentum4V1),0.0,0.0,0.0)
	if (tBoost1.get_momentum_in_CMF() == tExpected):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing get_momentum_in_CMF."
	assertions.pause(__name__)

	##Test forward and reverse boosts of CM to CMF:##
	print "\n--------------------------------------------------\n"
	print "Testing forward and reverse boosts of CM momentum vector to CMF:\n"
	print "CMF momentum vector = " , tMomentum4V1
	print "Momentum in = " , tMomentum4V1
	print "boost() = " , tBoost1.boost(tMomentum4V1)
	print "inverse_boost(boost()) = " , tBoost1.inverse_boost(tBoost1.boost(tMomentum4V1))
	if (tBoost1.inverse_boost(tBoost1.boost(tMomentum4V1)) == tMomentum4V1):
		print "\nReturned the original fourVector: Test successful!\n"
	else:
		print "\nTest Failed\n"
	print "boost(boost()) = " , tBoost1 * (tBoost1 * tMomentum4V1)
	print "inverse_boost(boost(boost())) = " , tBoost1 / (tBoost1 * (tBoost1 * tMomentum4V1))
	print "inverse_boost(inverse_boost(boost(boost()))) = " , tBoost1 / (tBoost1 / (tBoost1 * (tBoost1 * tMomentum4V1)))
	if ( (tBoost1 / (tBoost1 / (tBoost1 * (tBoost1 * tMomentum4V1)))) == tMomentum4V1):
		print "\nReturned the original fourVector: Test successful!"
	else:
		print "\nTest Failed"
	print "\nFinished testing forward and reverse boosts of CM momentum vector to CMF."
	assertions.pause(__name__)

	##Test forward and reverse boosts of arbitary (E,p) to CMF:##
	print "\n--------------------------------------------------\n"
	print "Testing forward and reverse boosts of arbitary (E,p) to CMF:\n"
	tCount = 0
	for tFourVector in [tFourVectors['a'],tFourVectors['b']]:
		print "~~~~~~~~\nBoosting" , tFourVector , "\ngives\n" , tBoost1 * tFourVector
		print "and then the inverse boost of that gives:\n" , tBoost1 / (tBoost1 * tFourVector)
		if (tBoost1.inverse_boost(tBoost1.boost(tFourVector)) == tFourVector):
			print "\nReturned the original fourVector: Test successful!"
		else:
			print "\nTest Failed"
		if tCount < 1:
			tCount += 1
			assertions.pause(__name__)
	print "\nFinished testing forward and reverse boosts of arbitary (E,p) to CMF."
	assertions.pause(__name__)

	##Test sucessive different (E,p) boosts:##
	print "\n--------------------------------------------------\n"
	print "Testing sucessive different (E,p) boosts:\n"
	print "Boosting" , tMomentum4V2 , "into the frame of:"
	print tMomentum4V1 , "gives:\n", tBoost1.boost(tMomentum4V2)
	tBoost2 = lorentzBoost(tBoost1.boost(tMomentum4V2))
	print "Then boosting that into its frame gives: " , tBoost2 * (tBoost1 * tMomentum4V2)
	if ((tBoost2 * (tBoost1 * tMomentum4V2)) == tBoost2.get_momentum_in_CMF()):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing sucessive different (E,p) boosts."
	assertions.pause(__name__)

	##Test boost of components and that Pi^2 is invariant:##
	print "\n--------------------------------------------------\n"
	print "Testing boost of components and that Pi^2 is invariant:\n"
	print "The summed centre of mass frame has momentum:\n", tSumMomentums
	for tFourVector in tFourVectors:
		print "\nThe Pi:" , tFourVectors[tFourVector] , "\nboosted to:" , tBoost3 * tFourVectors[tFourVector]
		print "\n(Pi)^2 before was" , tMomentumSquares[tFourVector] , "and after was" , tBoostedMomentumSquares[tFourVector]
		if precision.check_numbers_equal(tMomentumSquares[tFourVector],tBoostedMomentumSquares[tFourVector]):
			print "\nThis is conserved: Test successful!"
		else:
			print "{0:.20f}".format(tMomentumSquares[tFourVector])
			print "{0:.20f}".format(tBoostedMomentumSquares[tFourVector])
			print "\nNot conserved: Test failed!"
			print "Possibly due to rounding errors?"
			#Should look into this at some point - Possible *10^-10 new precision limit?
		print "\n---------"
		assertions.pause(__name__)
	print "The sum of the boosted Pi is:\n" , tBoostedSumMomentums , "\nand we tExpected:\n" , tBoost3.get_momentum_in_CMF()
	tSumResult = tBoost3.get_momentum_in_CMF() * tBoost3.get_momentum_in_CMF()
	print "Sum(Pi)^2 before was:" , tBoostedSumMomentums*tBoostedSumMomentums , "and now is:" , tSumResult
	t1 = (tBoostedSumMomentums == tBoost3.get_momentum_in_CMF())
	t2 = precision.check_numbers_equal(tBoostedSumMomentums*tBoostedSumMomentums,tSumResult)
	if (t1 and t2):
		print "\nThis is conserved: Test successful!"
	else:
		print "\nNot conserved: Test failed!"
		print "Possibly due to rounding errors?"
		print "\n{0:.20f}".format(tBoostedSumMomentums*tBoostedSumMomentums)
		print "{0:.20f}".format(tSumResult)
		#Same error as above in this testing for Pi^2.
	print "\nFinished testing boost of components and that Pi^2 is invariant."
	assertions.pause(__name__)

	##Test using massless particles:##
	print "\n--------------------------------------------------\n"
	print "Testing using massless particles:\n"
	print "Create a boost using:\n", tMasslessCOM
	print "\nBoost:\n", tMassless4V1, "\nof mass:", tMassless4V1*tMassless4V1,". Gives:\n" , tMasslessBoost1*tMassless4V1
	print "\nBoosting back gives:\n", tMasslessBoost1 / (tMasslessBoost1*tMassless4V1)
	print "\n----------\n"
	print "\nBoost:\n", tMassless4V2, "\nof mass:", tMassless4V2*tMassless4V2,". Gives:\n" , tMasslessBoost1*tMassless4V2
	print "\nBoosting back gives:\n", tMasslessBoost1 / (tMasslessBoost1*tMassless4V2)
	assertions.pause(__name__)
	tSuccessful = True
	print "\n----------\n"
	for ti in range(3):
		if not precision.check_numbers_equal((tMasslessBoost1*tMassless4V1)[ti+1],(-1.0 * (tMasslessBoost1*tMassless4V2)[ti+1])):
			tSuccessful = False
			print "Test failed as momentums not equal and opposite in COM frame."
		elif not precision.check_numbers_equal((tMasslessBoost1 / (tMasslessBoost1*tMassless4V1))[ti+1],tMassless4V1[ti+1]):
			tSuccessful = False
			print "Test failed as momentums not equal and opposite in COM frame."
		elif not precision.check_numbers_equal((tMasslessBoost1 / (tMasslessBoost1*tMassless4V2))[ti+1],tMassless4V2[ti+1]):
			tSuccessful = False
			print "Test failed as momentums not equal and opposite in COM frame."
	if tSuccessful:
		print "Momentums oppose in COM frame and boosting back works..."
	print "\nChecking light-like before returns:" , check_lightlike(tMassless4V1) , "and:" , check_lightlike(tMassless4V2)
	print "and afterwards:" , check_lightlike(tMasslessBoost1*tMassless4V1) , "and:" , check_lightlike(tMasslessBoost1*tMassless4V2)
	if (check_lightlike(tMassless4V1) == check_lightlike(tMasslessBoost1*tMassless4V1)) and tSuccessful:
		if (check_lightlike(tMassless4V2) == check_lightlike(tMasslessBoost1*tMassless4V2)):
			print "\nLight-like nature preserved: Test successful!\n"
		else:
			print "\nLight-like nature NOT preserved: Test failed!\n"
	else:
		print "\nLight-like nature NOT preserved: Test failed!\n"
	print "\nFinished testing using massless particles."
	assertions.pause(__name__)

	##Test lorentzRotation class:##
	print "\n--------------------------------------------------\n"
	print "//////////////////////////////"
	print "Testing lorentzRotation class:"
	print "//////////////////////////////"
	assertions.pause(__name__)

	##__init__() and calculate_angles tested implicitly through the other functions.

	##Test lorentzRotation forward and backwards on initalising vector:##
	print "\n--------------------------------------------------\n"
	print "Testing lorentzRotation forward and backwards on initalising vector:\n"
	print "Initialising vector = " , tMomentum4V1
	print "Vector in = " , tMomentum4V1
	print "rotate() = " , tRotate1.rotate(tMomentum4V1)
	print "inverse_rotate(rotate()) = " , tRotate1.inverse_rotate(tRotate1.rotate(tMomentum4V1))
	if (tRotate1.inverse_rotate(tRotate1.rotate(tMomentum4V1)) == tMomentum4V1):
		if ((tRotate1.rotate(tMomentum4V1)[1] == 0.0) and (tRotate1.rotate(tMomentum4V1)[2] == 0.0)):
			print "\nReturned the original fourVector: Test successful!\n"
		else:
			print "\nBoosted result not only in z-direction: Test Failed!\n"
	else:
		print "\nDidn't return original result: Test Failed!"
	print "rotate(rotate()) = " , tRotate1 * (tRotate1 * tMomentum4V1)
	print "inverse_rotate(rotate(rotate())) = " , tRotate1 / (tRotate1 * (tRotate1 * tMomentum4V1))
	print "inverse_rotate(inverse_rotate(rotate(rotate()))) = " , tRotate1 / (tRotate1 / (tRotate1 * (tRotate1 * tMomentum4V1)))
	if (tRotate1 / (tRotate1 / (tRotate1 * (tRotate1 * tMomentum4V1))) == tMomentum4V1):
		print "\nReturned the original fourVector: Test successful!\n"
	else:
		print "\nDidn't return original result: Test Failed!\n"
	print "\nFinished testing lorentzRotation forward and backwards on initalising vector."
	assertions.pause(__name__)

	##Test lorentzRotation forward and backwards on random vectors:##
	print "\n--------------------------------------------------\n"
	print "Testing lorentzRotation forward and backwards on random vectors:\n"
	print "Initialising vector = " , tMomentum4V1
	for tFourVector in [tFourVectors['a'],tFourVectors['b']]:
		print "Rotating" , tFourVector , "\ngives\n" , tRotate1 * tFourVector
		print "and then the inverse rotation of that gives:\n" , tRotate1 / (tRotate1 * tFourVector)
		if (tRotate1.inverse_rotate(tRotate1.rotate(tFourVector)) == tFourVector):
			print "\nReturned the original fourVector: Test successful!\n"
		else:
			print "\nTest Failed!\n"
	print "\nFinished testing lorentzRotation forward and backwards on random vectors."
	assertions.pause(__name__)

	##Test lorentzRotation forward and backwards on Pi and -Pi:##
	print "\n--------------------------------------------------\n"
	print "Testing lorentzRotation forward and backwards on Pi and -Pi:\n"
	print "Initialising vector = " , tMomentum4V1
	print "Pi in = " , tMomentum4V1
	print "-Pi in = " , tMomentum4V1Reversed
	print "\nrotate(Pi) = " , tRotate1.rotate(tMomentum4V1)
	print "rotate(-Pi) = " , tRotate1.rotate(tMomentum4V1Reversed)
	print "\ninverse_rotate(rotate(Pi)) = " , tRotate1.inverse_rotate(tRotate1.rotate(tMomentum4V1))
	print "inverse_rotate(rotate(-Pi)) = " , tRotate1.inverse_rotate(tRotate1.rotate(tMomentum4V1Reversed))
	tExpected = tRotate1.rotate(tMomentum4V1).copy()
	for ti in range(3):
		tExpected[ti+1] *= -1.0
	if ((tRotate1 / (tRotate1 * tMomentum4V1)) == tMomentum4V1):
		if ((tRotate1 / (tRotate1 * tMomentum4V1Reversed)) == tMomentum4V1Reversed):
			if (tRotate1 * tMomentum4V1Reversed) == tExpected:
				print "\nTests successful!\n"
			else:
				print "\nNot opposite after rotating: Test Failed!\n"
		else:
			print "\nReversed not the same after rotating back: Test Failed!\n"
	else:
		print "\nFourvector not the same after rotating back: Test Failed!\n"
	print "\nFinished testing lorentzRotation forward and backwards on Pi and -Pi."
	assertions.pause(__name__)

	##Test lorentzBoost and lorentzRotation classes together:##
	print "\n--------------------------------------------------\n"
	print "//////////////////////////////////////////////////////////"
	print "Testing lorentzBoost and lorentzRotation classes together:"
	print "//////////////////////////////////////////////////////////"
	assertions.pause(__name__)

	##Test lorentzBoost and lorentzRotation classes together on Pi and Pj:##
	print "\n--------------------------------------------------\n"
	print "Testing lorentzBoost and lorentzRotation forward and backwards on Pi and Pj:\n"
	print "Vector one is:" , tFourVectors['a']
	print "Vector two is:" , tFourVectors['b']
	print "\nCOM frame used for boost is =" , tTogetherCOM
	print "\nboost(V1) =", tBoost4*tFourVectors['a']
	print "boost(V2) =", tBoost4*tFourVectors['b']
	print "\nRotate set up from boost(V1) as chosen 'to decay'."
	print "rotate(boost(V1)) =", tRotate2*(tBoost4*tFourVectors['a'])
	print "rotate(boost(V2)) =", tRotate2*(tBoost4*tFourVectors['b'])
	tRV1, tRV2 = tRotate2*(tBoost4*tFourVectors['a']) , tRotate2*(tBoost4*tFourVectors['b'])
	t1 = precision.check_numbers_equal(tRV1[3],-1.0* tRV2[3])
	t2 = precision.check_numbers_equal(tRV1[1],0.0)
	t3 = precision.check_numbers_equal(tRV1[2],0.0)
	t4 = precision.check_numbers_equal(tRV2[1],0.0)
	t5 = precision.check_numbers_equal(tRV2[2],0.0)
	if (t1 and t2 and t3 and t4 and t5):
		print "\nTest successful!\n"
	else:
		print "\nTest Failed!\n"
	print "\nFinished testing lorentzBoost and lorentzRotation forward and backwards on Pi and Pj."
	assertions.pause(__name__)

	##Test lorentzBoost and lorentzRotation classes together on Pi and Pj of same mass:##
	print "\n--------------------------------------------------\n"
	print "Testing lorentzBoost and lorentzRotation forward and backwards on Pi and Pj of same mass:\n"
	print "Vector one is:" , tTogether2V1
	print "Vector two is:" , tTogether2V2C
	print "Their masses are:", tTogether2V1*tTogether2V1, "and:" , tTogether2V2C*tTogether2V2C
	print "\nCOM frame used for boost is =\n" , tTogether2COM
	print "boost(V1) =\n", tBoost5*tTogether2V1
	print "boost(V2) =\n", tBoost5*tTogether2V2C
	print "Rotate set up from boost(V1) as chosen 'to decay'."
	print "rotate(boost(V1)) =\n", tRotate3*(tBoost5*tTogether2V1)
	print "rotate(boost(V2)) =\n", tRotate3*(tBoost5*tTogether2V2C)
	tRV1, tRV2 = tRotate3*(tBoost5*tTogether2V1) , tRotate3*(tBoost5*tTogether2V2C)
	if (c_n_e(tRV1[0],tRV2[0]) and c_n_e(tRV1[3],-1.0* tRV2[3]) and c_n_e(tRV1[1],0.0) and c_n_e(tRV1[2],0.0) and c_n_e(tRV2[1],0.0) and c_n_e(tRV2[2],0.0)):
		print "\nTest successful!\n"
	else:
		print "\nTest Failed!"
		print "\nPossibly just a precision error?"
	assertions.pause(__name__)
	print "inverse_rotate(rotate(boost(V1))) =\n", tRotate3/(tRotate3*(tBoost5*tTogether2V1))
	print "inverse_rotate(rotate(boost(V2))) =\n", tRotate3/(tRotate3*(tBoost5*tTogether2V2C))
	print "\ninverse_boost(inverse_rotate(rotate(boost(V1)))) =\n", tBoost5/(tRotate3/(tRotate3*(tBoost5*tTogether2V1)))
	print "inverse_boost(inverse_rotate(rotate(boost(V2)))) =\n", tBoost5/(tRotate3/(tRotate3*(tBoost5*tTogether2V2C)))
	tRV1, tRV2 = tBoost5/(tRotate3/(tRotate3*(tBoost5*tTogether2V1))) , tBoost5/(tRotate3/(tRotate3*(tBoost5*tTogether2V2C)))
	if ((tRV1 == tTogether2V1) and (tRV2 == tTogether2V2C)):
		print "\nTest successful!\n"
	else:
		print "\n", tTogether2V1
		print tTogether2V2C
		print tRV1
		print tRV2
		print "\nTest Failed!\n"
	print "\nFinished testing lorentzBoost and lorentzRotation forward and backwards on Pi and Pj of same mass."
	assertions.pause(__name__)

	##Test lorentzBoost and lorentzRotation classes together on massless Pi and Pj:##
	print "\n--------------------------------------------------\n"
	print "Testing lorentzBoost and lorentzRotation forward and backwards on massless Pi and Pj:\n"
	print "Vector one is:" , tTogetherRV1Massless
	print "Vector two is:" , tTogetherRV2Massless
	print "Their masses are:", tTogetherRV1Massless*tTogetherRV1Massless, "and:" , tTogetherRV2Massless*tTogetherRV2Massless
	print "\nCOM frame used for boost is =\n" , tTogetherCOMMassless
	print "boost(V1) =\n", tBoost6*tTogetherRV1Massless
	print "boost(V2) =\n", tBoost6*tTogetherRV2Massless
	print "Rotate set up from boost(V1) as chosen 'to decay'."
	print "rotate(boost(V1)) =\n", tRotate4*(tBoost6*tTogetherRV1Massless)
	print "rotate(boost(V2)) =\n", tRotate4*(tBoost6*tTogetherRV2Massless)
	tRV1, tRV2 = tRotate4*(tBoost6*tTogetherRV1Massless) , tRotate4*(tBoost6*tTogetherRV2Massless)
	tCondition1, tCondition2 = precision.check_numbers_equal(tRV1[0],tRV2[0]), precision.check_numbers_equal(tRV1[3],-1.0* tRV2[3])
	if (tCondition1 and tCondition2 and c_n_e(tRV1[1],0.0) and c_n_e(tRV1[2],0.0) and c_n_e(tRV2[1],0.0) and c_n_e(tRV2[2],0.0)):
		if (c_n_e(tRV1*tRV1,0.0) and c_n_e(tRV2*tRV2,0.0)):
			print "\nTest successful!\n"
		else:
			print"\nResult no longer time-like: Test failed!\n"
	else:
		print "\n~~~~~~~~~\n", tRV1
		print tRV2
		print "\nResulting vectors not correct: Test Failed!\n"
		print "Possibly just a precision error?"
	assertions.pause(__name__)
	print "inverse_rotate(rotate(boost(V1))) =\n", tRotate4/(tRotate4*(tBoost6*tTogetherRV1Massless))
	print "inverse_rotate(rotate(boost(V2))) =\n", tRotate4/(tRotate4*(tBoost6*tTogetherRV2Massless))
	print "\ninverse_boost(inverse_rotate(rotate(boost(V1)))) =\n", tBoost6/(tRotate4/(tRotate4*(tBoost6*tTogetherRV1Massless)))
	print "inverse_boost(inverse_rotate(rotate(boost(V2)))) =\n", tBoost6/(tRotate4/(tRotate4*(tBoost6*tTogetherRV2Massless)))
	tRV1 = tBoost6/(tRotate4/(tRotate4*(tBoost6*tTogetherRV1Massless)))
	tRV2 = tBoost6/(tRotate4/(tRotate4*(tBoost6*tTogetherRV2Massless)))
	if ((tRV1 == tTogetherRV1Massless) and (tRV2 == tTogetherRV2Massless)):
		print "\nTest successful!"
	else:
		print "\n~~~~~~~~~\n", tTogetherRV1Massless
		print tTogetherRV2Massless
		print tRV1
		print tRV2
		print "\nTest Failed!"
		print "Possibly just a precision error?"
	print "\nFinished testing lorentzBoost and lorentzRotation forward and backwards on massless Pi and Pj."
	assertions.pause(__name__)

	##Test lorentzBoost class:##
	print "\n--------------------------------------------------\n"
	print "/////////////////////////////"
	print "Testing BoostAndRotate class:"
	print "/////////////////////////////"
	assertions.pause(__name__)

	##boost_rotate and inverse_boost_rotate tested implictly via the * and / below.

	##Test boostAndRotation classe on massless Pi and Pj:##
	print "\n--------------------------------------------------\n"
	print "Testing boostAndRotation forward and backwards on massless Pi and Pj:\n"
	print "Vector one is:" , tTogetherRV1Massless
	print "Vector two is:" , tTogetherRV2Massless
	print "Their masses are:", tTogetherRV1Massless*tTogetherRV1Massless, "and:" , tTogetherRV2Massless*tTogetherRV2Massless
	print "\nCOM frame used for boost is =\n" , tTogetherCOMMassless
	print "Initiated with index 0 to recoil:"
	print "boost_rotate(V1) =\n", tBoostAndRotate0*tTogetherRV1Massless
	print "boost_rotate(V2) =\n", tBoostAndRotate0*tTogetherRV2Massless
	assertions.pause(__name__)
	print "\nInitiated with index 1 to recoil:"
	print "boost_rotate(V1) =\n", tBoostAndRotate1*tTogetherRV1Massless
	print "boost_rotate(V2) =\n", tBoostAndRotate1*tTogetherRV2Massless
	tRV10, tRV20 = tBoostAndRotate0*tTogetherRV1Massless , tBoostAndRotate0*tTogetherRV2Massless
	tRV11, tRV21 = tBoostAndRotate1*tTogetherRV1Massless , tBoostAndRotate1*tTogetherRV2Massless
	tCondition1, tCondition2 = precision.check_numbers_equal(tRV10[0],tRV20[0]), precision.check_numbers_equal(tRV10[3],-1.0* tRV20[3])
	tCondition3, tCondition4 = precision.check_numbers_equal(tRV11[0],tRV21[0]), precision.check_numbers_equal(tRV11[3],-1.0* tRV21[3])
	t1 = (tCondition1 and tCondition2 and c_n_e(tRV10[1],0.0) and c_n_e(tRV10[2],0.0) and c_n_e(tRV20[1],0.0) and c_n_e(tRV20[2],0.0))
	t2 = (tCondition3 and tCondition4 and c_n_e(tRV11[1],0.0) and c_n_e(tRV11[2],0.0) and c_n_e(tRV21[1],0.0) and c_n_e(tRV21[2],0.0))
	if (t1 and t2):
		if (c_n_e(tRV10*tRV10,0.0) and c_n_e(tRV20*tRV20,0.0) and c_n_e(tRV11*tRV11,0.0) and c_n_e(tRV21*tRV21,0.0)):
			print "\nTest successful!\n"
		else:
			print"\nResult no longer time-like: Test failed!\n"
	else:
		print "\n~~~~~~~~~\n", tRV10
		print tRV20
		print tRV11
		print tRV21
		print "\nResulting vectors not correct: Test Failed!\n"
		print "Possibly just a precision error?"
	assertions.pause(__name__)
	print "\nInitiated with index 0 to recoil:"
	print "inverse_boost_rotate(boost_rotate(V1)) =\n", tBoostAndRotate0/tRV10
	print "inverse_boost_rotate(boost_rotate(V2)) =\n", tBoostAndRotate0/tRV20
	assertions.pause(__name__)
	print "\nInitiated with index 1 to recoil:"
	print "inverse_boost_rotate(boost_rotate(V1)) =\n", tBoostAndRotate1/tRV11
	print "inverse_boost_rotate(boost_rotate(V2)) =\n", tBoostAndRotate1/tRV21
	tEV10, tEV20, tEV11, tEV21 = tBoostAndRotate0/tRV10, tBoostAndRotate0/tRV20, tBoostAndRotate1/tRV11, tBoostAndRotate1/tRV21
	t1 = ((tEV10 == tTogetherRV1Massless) and (tEV20 == tTogetherRV2Massless))
	t2 = ((tEV11 == tTogetherRV1Massless) and (tEV21 == tTogetherRV2Massless))
	if (t1 and t2):
		print "\nTest successful!"
	else:
		print "\n~~~~~~~~~\n", tTogetherRV1Massless
		print tTogetherRV2Massless
		print tEV10
		print tEV20
		print tEV11
		print tEV21
		print "\nTest Failed!"
		print "Possibly just a precision error?"
	print "\nFinished testing boostAndRotation forward and backwards on massless Pi and Pj."
	assertions.pause(__name__)

	##Done testing:##
	print "\n---------------------------------------------\n"
	print "/////////////////////////////////"
	print "Finished checking lorentz module!"
	print "/////////////////////////////////"