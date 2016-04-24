####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module for handling dipole showering kinematics."""

##Import required modules:##
import math
import random
import assertions
import precision
import counters
import dataLoggers
import fourVectors
import lorentz

print "\n//////////////////////////"
print "Loading kinematics module:"
print "//////////////////////////\n"

##Set up data loggers:##
e1s = dataLoggers.dataLogger()
e2s = dataLoggers.dataLogger()
e3s = dataLoggers.dataLogger()
tX1s = dataLoggers.dataLogger()
tX3s = dataLoggers.dataLogger()

##Functions:##
	
def Sijk(listOfMomentumFourVectors):
	"""A function to calculate the invariant quantity Sijk for Pi's in the given list."""
	assert type(listOfMomentumFourVectors) == list
	__sqrtSijk = fourVectors.fourVector(0.0,0.0,0.0,0.0)
	for __Pi in listOfMomentumFourVectors:
		assert fourVectors.check_is_fourVector(__Pi)
		assert __Pi.__nonzero__()
		__sqrtSijk += __Pi
	__Sikj = __sqrtSijk * __sqrtSijk
	return __Sikj

def Pperp_squared(p1,p2,p3): ##Using p1' + p3' -> p1 + p2 + p3
	"""A function to calculate the lorentz invariant transverse momentum squared given the correct P1,2,3 vectors."""
	assert (fourVectors.check_is_fourVector(p1) and fourVectors.check_is_fourVector(p2) and fourVectors.check_is_fourVector(p3))
	assert (p1.__nonzero__() and p2.__nonzero__() and p3.__nonzero__())
	__S12, __S23, __S123 = Sijk([p1,p2]), Sijk([p2,p3]), Sijk([p1,p2,p3])
	__pPerpSquared = __S12 * __S23 / __S123
	return __pPerpSquared

def rapidity(p1,p2,p3): ##Using p1' + p3' -> p1 + p2 + p3
	"""A function to calculate the rapidity given the correct P1,2,3 vectors."""
	assert (fourVectors.check_is_fourVector(p1) and fourVectors.check_is_fourVector(p2) and fourVectors.check_is_fourVector(p3))
	assert (p1.__nonzero__() and p2.__nonzero__() and p3.__nonzero__())
	__S12, __S23 = Sijk([p1,p2]), Sijk([p2,p3])
	__y = 0.5*math.log(__S23/__S12)
	return __y

def E_to_S123(E):
	"""A function to return S123 from E."""
	assert assertions.all_are_numbers([E])
	return 4.0*E*E

def S123_to_E(S123):
	"""A function to return E from S123."""
	assert assertions.all_are_numbers([S123])
	return math.sqrt(S123)/2.0

##~~~~~~~~~~~~~~~~~~~##

def get_random_phi():
	"""A function to return a uniformly random value of phi in the range 0 to 2Pi."""
	__randomNumber = random.random()
	__randomPhi = 2.0 * math.pi * __randomNumber
	return __randomPhi

def get_random_theta():
	"""A function to return a uniformly random value of theta in the range 0 to Pi."""
	__randomNumber = random.random()
	__randomPhi = math.pi * __randomNumber
	return __randomPhi

def calculate_E1(S123,Pperp,y):
	"""A function to calculate a value of E1 following a dipole splitting."""
	assert assertions.all_are_numbers([S123,Pperp,y])
	__E1 = 0.5*(math.sqrt(S123) - (Pperp*math.exp(y)))
	return __E1

def calculate_E2(S123,Pperp,y):
	"""A function to calculate a value of E2 following a dipole splitting."""
	assert assertions.all_are_numbers([S123,Pperp,y])
	__E2 = Pperp*math.cosh(y)
	return __E2

def calculate_E3(S123,Pperp,y):
	"""A function to calculate a value of E3 following a dipole splitting."""
	assert assertions.all_are_numbers([S123,Pperp,y])
	__E3 = 0.5*(math.sqrt(S123) - (Pperp*math.exp((-1.0)*y)))
	return __E3

def calculate_kPerp(S123,Pperp,y):
	"""A function to calculate the actual transverse momentum following a dipole splitting."""
	assert assertions.all_are_numbers([S123,Pperp,y])
	__E1, __E2, __E3 = calculate_E1(S123,Pperp,y), calculate_E2(S123,Pperp,y), calculate_E3(S123,Pperp,y)
	assert (__E2 < (__E1 + __E3)) ##+ve y -> E3>E1 but -ve y -> E1>E3.
	if ((__E2 > __E1) or (__E2 > __E3)): ##Too many of these would be a problem so track occurance.
		counters.kPerpProdWarningCounter.count()
	__term2 = (((__E1*__E1) - (__E2*__E2) + (__E3*__E3))/(2.0*__E3))
	kPerpSquared = (__E1*__E1) - (__term2*__term2)
	return math.sqrt(kPerpSquared)

def calculate_x1(S123,Pperp,y):
	"""A function to return the momentum fraction X1 (e.g. Xq)."""
	assert assertions.all_are_numbers([S123,Pperp,y])
	__X1 = 1.0 - (Pperp*math.exp(y)/math.sqrt(S123))
	return __X1

def calculate_x2(S123,Pperp,y):
	"""A function to return the momentum fraction X2 (e.g. Xg or Xp)."""
	assert assertions.all_are_numbers([S123,Pperp,y])
	__X2 = (2.0*Pperp*math.cosh(y))/math.sqrt(S123)
	return __X2

def calculate_x3(S123,Pperp,y):
	"""A function to return the momentum fraction X3 (e.g. Xq_bar)."""
	assert assertions.all_are_numbers([S123,Pperp,y])
	__X3 = 1.0 - (Pperp*math.exp((-1.0)*y)/math.sqrt(S123))
	return __X3

def calculate_split_ps(S123,Pperp,y):
	"""A function to calculate the new p1,p2,p3 four-vectors following a dipole splitting."""
	##Letting p2 take +Pperp and p1 take -Pperp by convention.
	assert assertions.all_are_numbers([S123,Pperp,y])
	assert ((S123 > 0.0) and (Pperp > 0.0))
	__E1, __E2 = calculate_E1(S123,Pperp,y), calculate_E2(S123,Pperp,y)
	__E3, __kPerp = calculate_E3(S123,Pperp,y), calculate_kPerp(S123,Pperp,y)
	##Store all energies produced for plotting graphs after showering:
	e1s.store(__E1)
	e2s.store(__E2)
	e3s.store(__E3)
	__randomPhi, __newP1 = get_random_phi(), fourVectors.fourVector(0.0,0.0,0.0,0.0)
	__newP2, __newP3 = fourVectors.fourVector(0.0,0.0,0.0,0.0), fourVectors.fourVector(0.0,0.0,0.0,0.0)
	assert (0.0 < __randomPhi and __randomPhi < 2.0 * math.pi)
	##Populate the first vector:
	__newP1[0] = __E1
	__newP1[1] = (-1.0) * __kPerp * math.cos(__randomPhi)
	__newP1[2] = (-1.0) * __kPerp * math.sin(__randomPhi)
	__newP1[3] = math.sqrt((__E1*__E1) - (__kPerp*__kPerp))
	assert ((__newP1[0] > 0.0) and (__newP1[3] > 0.0))
	##Populate the second vector:
	__newP2[0] = __E2
	__newP2[1] = __kPerp * math.cos(__randomPhi)
	__newP2[2] = __kPerp * math.sin(__randomPhi)
	__newP2[3] = math.sqrt((__E2*__E2) - (__kPerp*__kPerp))
	assert (__newP2[0] > 0.0)
	##Populate the third vector. Entries 1 & 2 already 0 as required:
	__newP3[0] = __E3
	__newP3[3] = (-1.0) *__E3
	assert ((__newP3[0] > 0.0) and (__newP3[3] < 0.0))
	##Log x's for a produced particle for plotting after showering:
	tX1s.store(calculate_x1(S123,Pperp,y))
	tX3s.store(calculate_x3(S123,Pperp,y))
	##Check the direction of the produced particle is consistent with parallel momentum conservation.
	##This allows for different recoil configurations.
	if precision.check_numbers_equal(0.0,__newP1[3]+__newP2[3] + __newP3[3]):
		return __newP1, __newP2, __newP3
	elif precision.check_numbers_equal(0.0,__newP1[3] - __newP2[3] + __newP3[3]):
		__newP2[3] = -__newP2[3]
		return __newP1, __newP2, __newP3
	elif precision.check_numbers_equal(0.0,-__newP1[3] + __newP2[3] + __newP3[3]):
		__newP1[3] = -__newP1[3]
		return __newP1, __newP2, __newP3

def fix_energy_difference(difference,v2):
	"""A function to fix any energy difference between the two particles in and three out caused by code errors."""
	assert (type(difference) == float)
	assert fourVectors.check_is_fourVector(v2)
	assert v2.__nonzero__()
	##Place all of the energy change required into the new particle's vector, v2.
	__magIn = v2.calculate_cartesian_magnitude()
	__pFactor = math.sqrt((difference/(__magIn*__magIn)) + 1.0)
	__nV2 = v2.copy()
	__nV2 *= __pFactor
	__nV2[0] = math.sqrt((v2[0]*v2[0]) + difference)
	return __nV2

def produce_new_vectors(v1,v3,PperpSquared,y):
	"""A function to boost two vectors, perform the kinematics and return the tResults after boosting back."""
	##Requires v3 to be the recoiling vector, consistent with all kinematics in this project.
	assert (fourVectors.check_is_fourVector(v1) and fourVectors.check_is_fourVector(v3))
	assert (v1.__nonzero__() and v3.__nonzero__())
	assert assertions.all_are_numbers([PperpSquared,y])
	__energyIn = v1[0] + v3[0] ##This is in the 'lab frame' as these are never actually boosted in the code.
	__S123, __Pperp = Sijk([v1.copy(),v3.copy()]), math.sqrt(PperpSquared)
	__boostRotate = lorentz.boostAndRotate(v1.copy(),v3.copy(),1) ##Given v3 as vector to recoil so need index 1 here.
	__nBRV1, __nBRV2, __nBRV3 = calculate_split_ps(__S123,__Pperp,y)
	assert precision.check_numbers_equal(0.0,__nBRV1[3]+__nBRV2[3]+__nBRV3[3])
	__nV1, __nV2, __nV3 = __boostRotate/__nBRV1, __boostRotate/__nBRV2, __boostRotate/__nBRV3
	##Must also check energy out in the 'lab' frame too as not a Lorentz invariant quantity.
	__energyOut = __nV1[0] + __nV2[0] + __nV3[0]
	__energyDifference = __energyIn - __energyOut
	__nV2 = fix_energy_difference(__energyDifference,__nV2)
	##Check conservation rules:
	assert precision.check_numbers_equal(__energyIn,__energyOut)
	assert precision.check_numbers_equal(__S123,Sijk([__nV1.copy(), __nV2.copy(), __nV3.copy()]))
	##No check for any momentums here as no longer have to sum to zero in this frame.
	return __nV1, __nV2, __nV3

##Classes:##

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##
	from matplotlib import pyplot
	import sudakovs

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "//////////////////////////"
	print "Testing kinematics module:"
	print "//////////////////////////"
	assertions.pause(__name__)
	

	##Setup here:##
	print "\nGenerating test values..."
	tS123 = 0.0 ##Needed to start while loop.
	while tS123 < 10.0:
		##Generate random test four-vectors:
		tFourVectors = {'a':None,'b':None,'c':None,'d':None}
		##Prevent the random vectors being the same as can't boost into frame of massless particle.
		tVector1, tVector2 = fourVectors.fourVector(0,0,0,0), fourVectors.fourVector(0,0,0,0)
		while ((tVector1 == tVector2) or (not fourVectors.check_different_direction(tVector1,tVector2))):
			for tFourVector in tFourVectors:
				##Random range set so as to produce a vector slower than the speed of light.
				tX1 = random.randrange(2,100)
				tX2 = random.randrange(2,100)
				tX3 = random.randrange(2,100)
				##Make them massless.
				tX0 = math.sqrt((tX1*tX1) + (tX2*tX2) +(tX3*tX3))
				tFourVectors[tFourVector] = fourVectors.fourVector(tX0,tX1,tX2,tX3)
			tExpectedSab = tFourVectors['a']*tFourVectors['a'] + tFourVectors['b']*tFourVectors['b'] 
			tExpectedSab += 2.0*(tFourVectors['a']*tFourVectors['b'])
			tExpectedScd = tFourVectors['c']*tFourVectors['c'] + tFourVectors['d']*tFourVectors['d'] 
			tExpectedScd += 2.0*(tFourVectors['c']*tFourVectors['d'])
			tVector1 = tFourVectors['a'].copy()
			tVector2 = tFourVectors['b'].copy()
		#~~~~~~#
		print "Produced test vectors:"
		print tVector1 , "mass:", tVector1*tVector1
		print tVector2 , "mass:", tVector2*tVector2
		#~~~~~~#
		tCOM = tVector1 + tVector2
		tBoost1 = lorentz.lorentzBoost(tCOM)
		tP1NotRotated = tBoost1 * tVector1
		tP2NotRotated = tBoost1 * tVector2
		##Take a to be q and b to be qBar
		tRotation1 = lorentz.lorentzRotation(tP1NotRotated)
		tP1Before = tRotation1 * tP1NotRotated
		tP3Before = tRotation1 * tP2NotRotated
		#~~~~~~#
		print "Produced boosted test vectors:"
		print tP1Before, "mass:", tP1Before*tP1Before
		print tP3Before, "mass:", tP3Before*tP3Before
		#~~~~~~#
		tS123 = Sijk([tP1Before,tP3Before])
	tE = S123_to_E(tS123)
	tPperpSquared, tY = None, None
	print "Using S123 =", tS123
	while (tPperpSquared == None):
		tPperpSquared, tY, tProcCode = sudakovs.solve(tS123,tS123,1,-1) ##QQbar -> g/p
	tPperp = math.sqrt(tPperpSquared)
	tcalcX1 = calculate_x1(tS123,tPperp,tY)
	tcalcX2 = calculate_x2(tS123,tPperp,tY)
	tcalcX3 = calculate_x3(tS123,tPperp,tY)
	print "Produced energy fractions:"
	print "X1:", tcalcX1
	print "X3:", tcalcX3
	tY2 = 0.0
	tPperp2 = 0.0
	tP1After, tP2After, tP3After = calculate_split_ps(tS123,tPperp,tY)
	tS123After = Sijk([tP1After,tP2After,tP3After])
	#~~~~~~#
	tExpectedE1 = 0.5*(math.sqrt(tS123) - (tPperp*math.exp(tY)))
	tExpectedE2 = tPperp*math.sinh(2.0*tY)/(2.0*math.sinh(tY))
	tExpectedE3 = 0.5*(math.sqrt(tS123) - (tPperp*math.exp((-1.0)*tY)))
	tExpectedKPerpSquared = (tExpectedE1**2.0) - ((((tExpectedE1**2.0) - (tExpectedE2**2.0) + (tExpectedE3**2.0))/(2.0*tExpectedE3))**2.0)
	tExpectedKPerp = math.sqrt(tExpectedKPerpSquared)
	#~~~~~~#
	tKPerp1 = -1.0* math.sqrt((tP1After[1]*tP1After[1]) + (tP1After[2]*tP1After[2]))
	tKPerp2 = math.sqrt((tP2After[1]*tP2After[1]) + (tP2After[2]*tP2After[2]))
	c_n_e = precision.check_numbers_equal
	#~~~~~~#
	tEDiffVector = fourVectors.fourVector(math.sqrt(1.2*1.2 + 4.5*4.5 + 3.2*3.2),1.2,4.5,3.2)
	tEdiff = 5.4
	tEDiffVecOut = fix_energy_difference(tEdiff,tEDiffVector)
	print "Generating test variables completed"
	assertions.pause(__name__)

	##Test Sijk function:##
	print "\n--------------------------------------------------\n"
	print "Testing Sijk function:\n"
	for tFourVector in tFourVectors:
		print "Momentum-" + tFourVector + " is:\n" , tFourVectors[tFourVector]
	assertions.pause(__name__)
	print "\nFor a & b, Sijk returns:" , Sijk([tFourVectors['a'],tFourVectors['b']])
	print "Compared to the expected value of:" , tExpectedSab
	if precision.check_numbers_equal(Sijk([tFourVectors['a'],tFourVectors['b']]),tExpectedSab):
		print "\nFirst test, equal within the code precision: Test successful!"
	else:
		print "\nFirst test, not equal within the code precision: Test failed!"
		print "Possibly just a precision error?\n"
		print "{0:.20f}".format(Sijk([tFourVectors['a'],tFourVectors['b']]))
		print "{0:.20f}".format(tExpectedSab)
	print "\nFor a & b, Sijk returns:" , Sijk([tFourVectors['c'],tFourVectors['d']])
	print "Compared to the expected value of:" , tExpectedScd
	if precision.check_numbers_equal(Sijk([tFourVectors['c'],tFourVectors['d']]),tExpectedScd):
		print "\nSecond test, equal within the code precision: Test successful!"
	else:
		print "\nSecond test, not equal within the code precision: Test failed!"
		print "Possibly just a precision error?\n"
		print "{0:.20f}".format(Sijk([tFourVectors['c'],tFourVectors['d']]))
		print "{0:.20f}".format(tExpectedScd)
	print "\nFor all four Sijk returns:" , Sijk([tFourVectors['a'],tFourVectors['b'],tFourVectors['c'],tFourVectors['d']])
	print "\nFinished testing Sijk function."
	assertions.pause(__name__)

	##Test Pperp_squared function:##
	print "\n--------------------------------------------------\n"
	print "Testing Pperp_squared function:\n"
	for tFourVector in tFourVectors:
		print "Momentum-" + tFourVector + " is:\n" , tFourVectors[tFourVector]
	tResult = Pperp_squared(tFourVectors['a'],tFourVectors['b'],tFourVectors['c'])
	print "\nUsing a,b,c, the function returns:" , tResult
	tShouldBe = Sijk([tFourVectors['a'],tFourVectors['b']]) * Sijk([tFourVectors['b'],tFourVectors['c']])
	tShouldBe /= Sijk([tFourVectors['a'],tFourVectors['b'],tFourVectors['c']])
	if precision.check_numbers_equal(tResult,tShouldBe):
		print "\nAs expected: Test successful!"
	else:
		print "\nTest unsuccessful!"
	print "\nFinished testing Pperp_squared function."
	assertions.pause(__name__)

	##Test rapidity function:##
	print "\n--------------------------------------------------\n"
	print "Testing rapidity function:\n"
	for tFourVector in tFourVectors:
		print "Momentum-" + tFourVector + " is:\n" , tFourVectors[tFourVector]
	tResult = rapidity(tFourVectors['a'],tFourVectors['b'],tFourVectors['c'])
	print "\nUsing a,b,c, the function returns:" , tResult
	print "Sbc =", Sijk([tFourVectors['b'],tFourVectors['c']])
	print "Sab =", Sijk([tFourVectors['a'],tFourVectors['b']])
	tShouldBe = 0.5*math.log(Sijk([tFourVectors['b'],tFourVectors['c']]) / Sijk([tFourVectors['a'],tFourVectors['b']]))
	if precision.check_numbers_equal(tResult,tShouldBe):
		print "\nAs expected: Test successful!"
	else:
		print "\nTest unsuccessful! Could be incorrect test vectors giving S12 or S23 as -ve?"
	print "\nFinished testing rapidity function."
	assertions.pause(__name__)

	##Test E_to_S123 and S123_to_E functions:##
	print "\n--------------------------------------------------\n"
	print "Testing E_to_S123 and S123_to_E functions:\n"
	print "Using E =", tE, "and S123 =" , tS123
	print "E_to_S123(E) returns:", E_to_S123(tE)
	print "S123_to_E(S123) returns:", S123_to_E(tS123)
	if (precision.check_numbers_equal(tS123,E_to_S123(tE)) and precision.check_numbers_equal(tE,S123_to_E(tS123))):
		print "\nAs expected: Test successful!"
	else:
		print "\nTest unsuccessful!"
	print "\nFinished testing E_to_S123 and S123_to_E functions."
	assertions.pause(__name__)

	##Test get_random_phi function:##
	print "\n--------------------------------------------------\n"
	print "Testing get_random_phi function:\n"
	tRangeMultiplier = [1,1000]
	tSums = [[0],[0]]
	tDifferences = []
	for tSumsIndex, m in enumerate(tRangeMultiplier):
		tResults = []
		tSumPhi = 0.0
		tiRange = 10000*m
		for ti in range(10000):
			tResults.append(get_random_phi())
			tSumPhi += tResults[-1]
			tSums[tSumsIndex].append(tSumPhi/(ti+1))
		tAveragePhi = tSumPhi / len(tResults)
		print "Using" , tiRange , "function calls:"
		print "The average of the tResults is:" , tAveragePhi
		print "Which is:" , tAveragePhi - math.pi , "from the expected, Pi."
		tDifferences.append(tAveragePhi - math.pi)
		pyplot.figure()
		pyplot.plot(tResults,color = "blue", label = r"$Current\ value$")
		pyplot.plot(tSums[tSumsIndex], linewidth = 2, color = 'black',label = r"$Average\ value$")
		pyplot.axhline(2.0*math.pi, linestyle = '--', color = 'red')
		pyplot.axhline(math.pi, linestyle = '--', color = 'red')
		pyplot.axhline(tAveragePhi, linestyle = '--', color = 'green')
		pyplot.title(r"$Random\ values\ of\ \phi\ for\ " + str(tiRange) + r"\ iterations$")
		pyplot.ylabel(r"$Random\ \phi\ (rad)$")
		pyplot.xlabel(r"$Iteration$")
		pyplot.legend()
		assertions.show_graph()
	if (tDifferences[1]**2 < tDifferences[0]**2):
		print "\nAverage closer to pi with more calls:: Test successful!"
	else:
		print "\nAverage convergence test unsuccessful, check on graphs!"
	print "\nFinished testing get_random_phi function."
	assertions.pause(__name__)

	##Test get_random_theta function:##
	print "\n--------------------------------------------------\n"
	print "Testing get_random_theta function:\n"
	tRangeMultiplier = [1,1000]
	tSums = [[0],[0]]
	tDifferences = []
	for tSumsIndex, tm in enumerate(tRangeMultiplier):
		tResults = []
		tSumTheta = 0.0
		tiRange = 10000*tm
		for ti in range(10000):
			tResults.append(get_random_theta())
			tSumTheta += tResults[-1]
			tSums[tSumsIndex].append(tSumTheta/(ti+1))
		tAverageTheta = tSumTheta / len(tResults)
		print "Using" , tiRange , "function calls:"
		print "The average of the tResults is:" , tAverageTheta
		print "Which is:" , tAverageTheta - math.pi/2.0 , "from the expected, Pi/2."
		tDifferences.append(tAverageTheta - math.pi/2.0)
		pyplot.figure()
		pyplot.plot(tResults,color = "blue", label = r"$Current\ value$")
		pyplot.plot(tSums[tSumsIndex], linewidth = 2, color = 'black',label = r"$Average\ value$")
		pyplot.axhline(math.pi, linestyle = '--', color = 'red')
		pyplot.axhline(math.pi/2.0, linestyle = '--', color = 'red')
		pyplot.axhline(tAverageTheta, linestyle = '--', color = 'green')
		pyplot.title(r"$Random\ values\ of\ \theta\ for\ " + str(tiRange) + r"\ iterations$")
		pyplot.ylabel(r"$Random\ \theta\ (rad)$")
		pyplot.xlabel(r"$Iteration$")
		pyplot.legend()
		assertions.show_graph()
	if (tDifferences[1]**2 < tDifferences[0]**2):
		print "\nAverage closer to pi/2 with more calls:: Test successful!"
	else:
		print "\nAverage convergence test unsuccessful, check on graphs!"
	print "\nFinished testing get_random_theta function."
	assertions.pause(__name__)

	##Test calculate_E1, calculate_E2 and calculate_E3 functions:##
	print "\n--------------------------------------------------\n"
	print "Testing calculate_E1, calculate_E2 and calculate_E3 functions:\n"
	print "Using S123 =" , tS123
	print "Pperp =" , tPperp
	print "and y =" , tY
	print "calculate_E1 returns:", calculate_E1(tS123,tPperp,tY)
	print "calculate_E2 returns:", calculate_E2(tS123,tPperp,tY)
	print "calculate_E3 returns:", calculate_E3(tS123,tPperp,tY)
	tExpected1 = 0.5*(math.sqrt(tS123) - tPperp*math.exp(tY))
	tExpected2 = tPperp*math.cosh(tY)
	tExpected3 = 0.5*(math.sqrt(tS123) - tPperp*math.exp(-tY))
	tCheck1 = c_n_e(tExpected1,calculate_E1(tS123,tPperp,tY))
	tCheck2 = c_n_e(tExpected2,calculate_E2(tS123,tPperp,tY))
	tCheck3 = c_n_e(tExpected3,calculate_E3(tS123,tPperp,tY))
	assertions.pause(__name__)
	print "\nUsing S123 =" , tS123
	print "Pperp =" , tPperp
	print "and y =" , tY2
	print "calculate_E1 returns:", calculate_E1(tS123,tPperp,tY2)
	print "calculate_E2 returns:", calculate_E2(tS123,tPperp,tY2)
	print "calculate_E3 returns:", calculate_E3(tS123,tPperp,tY2)
	tExpected4 = 0.5*(math.sqrt(tS123) - tPperp)
	tExpected5 = tPperp
	tExpected6 = 0.5*(math.sqrt(tS123) - tPperp)
	tCheck4 = c_n_e(tExpected4,calculate_E1(tS123,tPperp,tY2))
	tCheck5 = c_n_e(tExpected5,calculate_E2(tS123,tPperp,tY2))
	tCheck6 = c_n_e(tExpected6,calculate_E3(tS123,tPperp,tY2))
	assertions.pause(__name__)
	print "\nUsing S123 =" , tS123
	print "Pperp =" , tPperp2
	print "and y =" , tY
	print "calculate_E1 returns:", calculate_E1(tS123,tPperp2,tY)
	print "calculate_E2 returns:", calculate_E2(tS123,tPperp2,tY)
	print "calculate_E3 returns:", calculate_E3(tS123,tPperp2,tY)
	tExpected7 = 0.5*math.sqrt(tS123)
	tExpected8 = 0.0
	tExpected9 = 0.5*math.sqrt(tS123)
	tCheck7 = c_n_e(tExpected7,calculate_E1(tS123,tPperp2,tY))
	tCheck8 = c_n_e(tExpected8,calculate_E2(tS123,tPperp2,tY))
	tCheck9 = c_n_e(tExpected9,calculate_E3(tS123,tPperp2,tY))
	if (tCheck1 and tCheck2 and tCheck3 and tCheck4 and tCheck5 and tCheck6 and tCheck7 and tCheck8 and tCheck9):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing calculate_E1, calculate_E2 and calculate_E3 functions."
	assertions.pause(__name__)

	##Test calculate_kPerp function:##
	print "\n--------------------------------------------------\n"
	print "Testing calculate_kPerp function:\n"
	print "Using S123 =" , tS123
	print "Pperp =" , tPperp
	print "and y =" , tY
	print "Calculate_kPerp returns:", calculate_kPerp(tS123,tPperp,tY)
	tE1, tE2 = calculate_E1(tS123,tPperp,tY), calculate_E2(tS123,tPperp,tY)
	tE3 = calculate_E3(tS123,tPperp,tY)
	term1 = (tE1*tE1-tE2*tE2+tE3*tE3)/(2.0*tE3)
	tExpected1 = math.sqrt(tE1*tE1 - term1*term1)
	print "And expected:", tExpected1
	if precision.check_numbers_equal(tExpected1,calculate_kPerp(tS123,tPperp,tY)):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing calculate_kPerp function."
	assertions.pause(__name__)

	##Test calculate_split_ps functions:##
	print "\n--------------------------------------------------\n"
	print "Testing calculate_split_ps functions:\n"
	print "Using S123 =", tS123 , "Pperp =", tPperp,  "y =", tY
	print "\nP1(before):" , tP1Before
	print "\nP3(before):" , tP3Before
	print "Their masses are:" , tP1Before*tP1Before , "and:" , tP3Before*tP3Before
	print "Calling calculate_split_ps..."
	tP1After , tP2After, tP3After = calculate_split_ps(tS123,tPperp,tY)
	print "It returned:"
	print "\nP1:", tP1After
	print "Which has mass:", tP1After * tP1After
	assertions.pause(__name__)
	print "\nP2:", tP2After
	print "Which has mass:", tP2After * tP2After
	assertions.pause(__name__)
	print "We expect kPerp:", tExpectedKPerp
	print "P1(after) gives:", tKPerp1
	print "Which implies phi was:", -1.0*math.atan2(tP1After[2],tPperp)
	print "P2(after) gives:", tKPerp2
	print "Which implies phi was:", math.atan2(tP2After[2],tPperp)
	if (precision.check_numbers_equal(-1.0*tKPerp1,tExpectedKPerp) and precision.check_numbers_equal(tKPerp2,tExpectedKPerp)):
		print "kPerp all equal within the code precision!"
	else:
		print "kPerp NOT equal within the code precision!"
		assert False
	if (precision.check_numbers_equal(tP1After[1],-1.0*tP2After[1]) and precision.check_numbers_equal(tP1After[2],-1.0*tP2After[2])):
		print "Individual x and y kPerps conserve momentum within the code precision!"
	else:
		print "Individual x and y kPerps DON'T conserve momentum within the code precision!"
		assert False
	assertions.pause(__name__)
	print "\nIt also returns P3:", tP3After
	print "Which has mass:", tP3After * tP3After
	if precision.check_numbers_equal(tP3After[0],-1.0*tP3After[3]):
		print "P3(after)[0] is -P3(after)[3] as expected: Test successful!"
	else:
		print "P3(after)[0] NOT -P3(after)[3] as expected: Test failed!"
		assert False
	tCheckE1 = precision.check_numbers_equal(tP1After[0],tExpectedE1)
	tCheckE2 = precision.check_numbers_equal(tP2After[0],tExpectedE2)
	tCheckE3 = precision.check_numbers_equal(tP3After[0],tExpectedE3)
	assertions.pause(__name__)
	if (tCheckE1 and tCheckE2 and tCheckE3):
		print "All energies as expected: Test successful!"
	else:
		print "All energies NOT as expected: Test failed!"
		assert False
	if precision.check_numbers_equal(tP1After[0] + tP2After[0] + tP3After[0],tP1Before[0] + tP3Before[0]):
		print "Energy conserved within the code precision!"
	else:
		print "Energy before:", tP1Before[0] + tP3Before[0]
		print "Energy after:", tP1After[0] + tP2After[0] + tP3After[0]
		print "Energy NOT conserved within the code precision!"
		assert False
	if precision.check_numbers_equal(tP1After[3] + tP2After[3] + tP3After[3],0.0):
		print "Pparallel conserved within the code precision!"
	else:
		print "Pparallel NOT conserved within the code precision!"
		assert False
	assertions.pause(__name__)
	print "S123 before was:" , tS123
	print "These give a new S123 of:" , tS123After
	if precision.check_numbers_equal(tS123,tS123After):
		print "These are equal: Test successful!"
	else:
		print "Test unsuccessful!"
		assert False
	assertions.pause(__name__)
	print "Used the values X1:", tcalcX1, "X3:", tcalcX3
	print "Calling calculate_x1, calculate_x2 and calculate_x3 gives:"
	print "X1:", calculate_x1(tS123,tPperp,tY), "X3:",  calculate_x3(tS123,tPperp,tY)
	print "Xi = 2Ei/Root(S123) gives X1:", tP1After[0]/tP1Before[0], "and tX3:",  tP3After[0]/tP3Before[0]
	print "Xi = PiParallel/PiParallelBefore gives X1:", tP1After[3]/tP1Before[3], "and tX3:",  tP3After[3]/tP3Before[3]
	tCheck1 = precision.check_numbers_equal(tcalcX1, calculate_x1(tS123,tPperp,tY))
	tCheck2 = precision.check_numbers_equal(tcalcX3, calculate_x3(tS123,tPperp,tY))
	tCheck3 = precision.check_numbers_equal(tcalcX1, tP1After[0]/tP1Before[0])
	tCheck4 = precision.check_numbers_equal(tcalcX3, tP3After[0]/tP3Before[0])
	tCheck5 = precision.check_numbers_equal(tcalcX1, tP1After[3]/tP1Before[3])
	tCheck6 = precision.check_numbers_equal(tcalcX3, tP3After[3]/tP3Before[3])
	if (tCheck1 and tCheck2):
		if (tCheck3 and tCheck4):
			print "Xq and Xq_bar all agree: Test successful!"
		else:
			print "First Xi definition used is NOT consistent: Test unsuccessful!"
			assert False
	else:
		print "Xq and Xq_bar function tResults do not match the X's used: Test unsuccessful!"
		assert False
	if (tCheck5 and tCheck6):
		print "\nSecond Xi definition is consistent. Unexpected!"
	else:
		print "\nSecond Xi definition is not consistent, as expected."
	print "\nFinished testing calculate_split_ps functions."
	assertions.pause(__name__)

	##Test calculate_x1, calculate_x2 and calculate_x3 functions:##
	print "\n--------------------------------------------------\n"
	print "Testing calculate_x1, calculate_x2 and calculate_x3 functions:\n"
	print "Using S123 =", tS123 , "Pperp =", tPperp,  "y =", tY
	print "X1 was:", tcalcX1, "and X3 was:", tcalcX3
	#testXq, testXq_bar = calculate_x1(tS123,tPperp,tY), calculate_x3(tS123,tPperp,tY)
	print "calculate_x1 returns:", tcalcX1
	print "calculate_x2 returns:", tcalcX2
	print "calculate_x3 returns:", tcalcX3
	print "They sum to:", tcalcX1 + tcalcX2 + tcalcX3
	print "Pperp squared is:", tPperp*tPperp
	print "Pperp squared = S123(1-Xq)(1-Xq_bar) is:", tS123*(1.0-tcalcX1)*(1.0-tcalcX3)
	print "Using the three tResulting vectors, Pperp squared is:", Pperp_squared(tP1After,tP2After,tP3After)
	if precision.check_numbers_equal(tPperp*tPperp,tS123*(1.0-tcalcX1)*(1.0-tcalcX3)):
		if precision.check_numbers_equal(tPperp*tPperp,Pperp_squared(tP1After,tP2After,tP3After)):
			print "\nAll as expected: Test successful!"
		else:
			print "\nError with Pperb^2 function vs Pperp^2: Test unsuccessful!"
			assert False
	else:
		print "\nError with S123(1-Xq)(1-Xq_bar) vs Pperp^2: Test unsuccessful!"
		assert False
	print "\nFinished testing calculate_x1, calculate_x2 and calculate_x3 functions."
	assertions.pause(__name__)

	##Test fix_energy_difference function:##
	print "\n--------------------------------------------------\n"
	print "Testing fix_energy_difference function:\n"
	print "Using", tEDiffVector
	print "which has mass", tEDiffVector*tEDiffVector, "and energy", tEDiffVector[0]
	print "Given an energy difference of", tEdiff
	print "returns:", tEDiffVecOut
	print "which has mass", tEDiffVecOut*tEDiffVecOut, "and energy", tEDiffVecOut[0]
	tCheck1 = precision.check_numbers_equal(tEDiffVecOut*tEDiffVecOut,0.0)
	tCheck2 = ((tEDiffVector != tEDiffVecOut) and (not fourVectors.check_different_direction(tEDiffVecOut,tEDiffVector)))
	if (tCheck1 and tCheck2):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
		assert False
		print "\nPossibly just a precision error?"
	print "\nFinished testing fix_energy_difference function."
	assertions.pause(__name__)

	##Test produce_new_vectors function:##
	print "\n--------------------------------------------------\n"
	print "Testing produce_new_vectors function:\n"
	print "Using Pperp =", tPperp,  "y =", tY
	print tP1Before
	print tP3Before
	print "\nCalling produce_new_vectors returned:"
	tNV1, tNV2, tNV3 = produce_new_vectors(tP1Before,tP3Before,tPperp,tY)
	print tNV1
	print tNV2
	print tNV3
	tM1, tM2, tM3 = tNV1*tNV1, tNV2*tNV2, tNV3*tNV3
	print "Which have masses:", tM1, tM2, tM3
	print "S123 before was:", Sijk([tP1Before,tP3Before])
	print "S123 after is:", Sijk([tNV1,tNV2,tNV3])
	tCheck1 = precision.check_numbers_equal(Sijk([tNV1,tNV2,tNV3]),Sijk([tP1Before,tP3Before]))
	tCheck2 = (precision.check_numbers_equal(tM1,0) and precision.check_numbers_equal(tM2,0) and precision.check_numbers_equal(tM3,0))
	if (tCheck1 and tCheck2):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
		assert False
		print "\nPossibly just a precision error?"
	print "\nFinished testing produce_new_vectors function."
	assertions.pause(__name__)
		
	##Done testing:##
	print "\n---------------------------------------------\n"
	print "////////////////////////////////////"
	print "Finished checking kinematics module!"
	print "////////////////////////////////////"