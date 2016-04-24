####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module for creating and handling four-vector objects."""

##Import required modules:##
import math
import assertions
import precision

print "\n///////////////////////////"
print "Loading fourVectors module:"
print "///////////////////////////\n"

##Functions:##

def check_is_fourVector(toCheck):
	"""A function to check for an instance of the fourVector class."""
	return isinstance(toCheck,fourVector)

def check_different_direction(toCheck1,toCheck2):
	"""A function to check whether two fourVectors point in different directions."""
	assert (check_is_fourVector(toCheck1) and check_is_fourVector(toCheck2))
	assert (toCheck1.__nonzero__() and toCheck2.__nonzero__())
	__uV1 = toCheck1.get_cartesian_unit_vector()
	__uV2 = toCheck2.get_cartesian_unit_vector()
	__zV = fourVector(0.0,0.0,0.0,0.0)
	if (__uV1 == __uV2):
		return False
	elif ((__uV1 == __zV) or (__uV2 == __zV)):
		return False
	else:
		return True

##Classes:##

class fourVector(object):
	"""A class for a four-vector."""
	
	def count_initially_empty(self):
		"""A function to count the number of empty inputs passed to an initiated four-vector object."""
		self.__numberInitiallyEmpty = 0
		for __i1, __tX1 in enumerate([self.__initialX0,self.__initialX1,self.__initialX2,self.__initialX3]):
			if __tX1 == None:
				self.__numberInitiallyEmpty += 1
				self.__positionsEmpty[__i1] = True

	def force_float_vector(self):
		"""A function to convert all vector integers to floats but leave doubles."""
		for __i2 in range(4):
			assert assertions.all_are_numbers([__i2])
			if (not assertions.check_float(self.__vector[__i2])):
				self.__vector[__i2] = assertions.force_float_number(self.__vector[__i2])

	def __init__(self, tX0=None, tX1=None, tX2=None, tX3=None):
		"""A function to initiate a four-vector object."""
		for __anX in [tX0, tX1, tX2, tX3]:
			assert ((__anX == None) or assertions.all_are_numbers([__anX]))
		self.__positionsEmpty = [False,False,False,False]
		self.__initialX0 , self.__initialX1 , self.__initialX2, self.__initialX3 = tX0 , tX1 , tX2 , tX3
		##Track un-filled four-vectors (using self.__isEmpty = T/F) to prevent empty vector appearing as zero vector.
		self.count_initially_empty()
		assert ((self.__numberInitiallyEmpty == 0) or (self.__numberInitiallyEmpty == 4))
		if (self.__numberInitiallyEmpty == 0):
			self.__isEmpty = False
			self.__initialX0, self.__initialX1 = self.__initialX0 , self.__initialX1
			self.__initialX2, self.__initialX3 = self.__initialX2 , self.__initialX3
			vector = [self.__initialX0,self.__initialX1,self.__initialX2,self.__initialX3]
			assert assertions.all_are_numbers(vector)
			self.__vector = vector
		elif self.__numberInitiallyEmpty == 4:
			self.__isEmpty = True
			self.__vector = [0.0, 0.0, 0.0, 0.0]
		self.force_float_vector()

	def copy(self):
		"""A function to create an exact copy of a four-vector."""
		__result1 = fourVector()
		for __i3 in range(4):
			__result1[__i3] = self.__vector[__i3]
		return __result1

	def update_empty(self):
		"""A function to check whether a four-vector is no longer empty."""
		if (all(__i4 == False for __i4 in self.__positionsEmpty)):
			self.__isEmpty = False

	def calculate_cartesian_magnitude(self):
		"""A function to return the magnitude of the three-vector in a four-vector."""
		assert self.__nonzero__()
		return math.sqrt((self[1] * self[1]) + (self[2] * self[2]) + (self[3] * self[3]))

	def get_cartesian_unit_vector(self):
		"""A function to return the cartesian unit three vector from a four-vector."""
		__magnitude = self.calculate_cartesian_magnitude()
		if (__magnitude == 0):
			return fourVector(0.0,0.0,0.0,0.0)
		__direction = fourVector(0.0,self[1],self[2],self[3])
		__direction /= __magnitude
		return __direction

	def __nonzero__(self):
		"""A function to determine whether the four-vector is a zero object (empty / not fully initiated)."""
		return not self.__isEmpty

	def __getitem__(self, i5):
		"""A function to get an entry in a four-vector object."""
		assert assertions.valid_four_vector_index(i5)
		assert self.__nonzero__()
		return self.__vector[i5]

	def __setitem__(self, i6, newValue):
		"""A function to assign an entry in a four-vector object."""
		assert assertions.valid_four_vector_index(i6)
		assert assertions.all_are_numbers([newValue])
		__newValue = assertions.force_float_number(newValue)
		self.__vector[i6] = __newValue
		self.__positionsEmpty[i6] = False
		self.update_empty()

	def __str__(self):
		"""A function to return a string of a four-vector object."""
		if (self.__isEmpty == True):
			return "[empty four-vector" + str(self.__positionsEmpty) + "]"
		else:
			return ("[" + str(self.__vector) + "]")

	def __repr__(self):
		"""A function to return a representation of a four-vector object."""
		if self.__isEmpty == True:
			return ("<four-vector([empty four-vector[" + str(self.__positionsEmpty) + "]])>")
		else:
			__vec = self.__vector
			return ("<four-vector([tX0: " + str(__vec[0]) + ", tX1: " + str(__vec[1]) + ", tX2: " + str(__vec[2]) + ", tX3: " + str(__vec[3]) + "])>")

	def __imul__(self,multiplier): ##Used for *=
		"""A function for multiplying a four-vector by a scalar."""
		assert assertions.all_are_numbers([multiplier])
		assert self.__nonzero__()
		__multiplier = assertions.force_float_number(multiplier)
		for __i7 in range(4):
			self[__i7] *= __multiplier
		return self

	def __mul__(fVA, fVB): ##Used for * to give Minkowski inner product.
		"""A function for calculating the Minkowski inner product of four-vector objects."""
		assert (check_is_fourVector(fVA) and check_is_fourVector(fVB))
		assert (fVA.__nonzero__() and fVB.__nonzero__())
		__result2 = ( fVA[0] * fVB[0] ) - ( (fVA[1] * fVB[1]) + (fVA[2] * fVB[2]) + (fVA[3] * fVB[3]) )
		return __result2

	def __idiv__(self,denominator): ##Used for /=
		"""A function for dividing a four-vector by a scalar."""
		assert assertions.all_are_numbers([denominator])
		assert assertions.safe_division(denominator)
		assert self.__nonzero__()
		__denominator = assertions.force_float_number(denominator)
		for __i8 in range(4):
			self[__i8] /= __denominator
		return self

	def __div__(self,denominator): ##Used for fVA / denominator.
		"""A function for dividing a four-vector by a scalar."""
		assert assertions.all_are_numbers([denominator])
		assert assertions.safe_division(denominator)
		assert self.__nonzero__()
		__denominator = assertions.force_float_number(denominator)
		__result3 = self.copy()
		__result3 /= __denominator
		return __result3

	def __iadd__(self,toAdd): ##Used for += on four-vectors.
		"""A function for adding another four-vector to itself."""
		assert check_is_fourVector(toAdd)
		assert (self.__nonzero__() and toAdd.__nonzero__())
		__result4 = fourVector()
		for __i9 in range(4):
			__result4[__i9] = self[__i9] + toAdd[__i9]
		return __result4

	def __add__(fVA, fVB): ##Used for a + b on four-vectors.
		"""A function for addition of four-vector objects."""
		assert (check_is_fourVector(fVA) and check_is_fourVector(fVB))
		assert (fVA.__nonzero__() and fVB.__nonzero__())
		__result5 = fVA.copy()
		__result5 += fVB
		return __result5

	def __isub__(self,toSubtract): ##Used for -= on four-vectors.
		"""A function for subtracting another four-vector from itself."""
		assert check_is_fourVector(toSubtract)
		assert (self.__nonzero__() and toSubtract.__nonzero__())
		__result6 = fourVector()
		for __i10 in range(4):
			__result6[__i10] = self[__i10] - toSubtract[__i10]
		return __result6

	def __sub__(fVA,fVB): ##Used for a - b on four-vectors.
		"""A function for subtracting one four-vector from another."""
		assert (check_is_fourVector(fVA) and check_is_fourVector(fVB))
		assert (fVA.__nonzero__() and fVB.__nonzero__())
		__result7 = fVA.copy()
		__result7 -= fVB
		return __result7

	def __abs__(self):
		"""A function to return the absolute value of a four-vector."""
		assert self.__nonzero__()
		return self * self

	def __eq__(self,other):
		"""A function to check whether two four-vectors are equal up to the code precision."""
		assert check_is_fourVector(other)
		assert (self.__nonzero__() and other.__nonzero__())
		for __i11 in range(4):
			if not precision.check_numbers_equal(self[__i11],other[__i11]):
				return False
		return True

	def __ne__(self,other):
		"""A function to check whether two four-vectors are not equal up to the code precision."""
		return not self.__eq__(other)

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##
	import random

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "///////////////////////////"
	print "Testing fourVectors module:"
	print "///////////////////////////"
	assertions.pause(__name__)
	
	##Setup here:##
	print "\nGenerating test values..."
	##Generate random test four-vectors:##
	tFourVectors = {'a':None,'b':None,'c':None,'d':None}
	for tFourVector in tFourVectors:
		tX0 = random.randrange(20)
		tX1 = random.randrange(20)
		tX2 = random.randrange(20)
		tX3 = random.randrange(20)
		tFourVectors[tFourVector] = fourVector(tX0,tX1,tX2,tX3)
	##~~~~~~~##
	tFloat = 2.30000
	tString = "string"
	##~~~~~~~##
	tSameDirection = tFourVectors['a'].copy()
	tSameDirection *= 2.5
	tOppositeDirection = tSameDirection.copy()
	tOppositeDirection *= -1.0
	tOneOppositeDirection = tSameDirection.copy()
	tOneOppositeDirection[2] = (-1.0) * tOneOppositeDirection[2]
	tOneZero = tSameDirection.copy()
	tOneZero[random.randrange(1,4)] = 0.0
	tTwoZero = tSameDirection.copy()
	tTwoZero[1], tTwoZero[2] = 0.0, 0.0
	tAllZero = tSameDirection.copy()
	tAllZero *= 0.0
	##~~~~~~~##
	tEmpty1 = fourVector()
	tEmpty1.count_initially_empty()
	tEmpty2 = fourVector(1,2,3,4)
	tEmpty2.count_initially_empty()
	##~~~~~~~##
	tMagnitude = fourVector(50,1.0,2.0,0.5)
	##~~~~~~~~#
	tUnitVector1 = fourVector(10,1.0,1.0,1.0)
	tUnitVector2 = fourVector(-5,2.0,-2.0,2.0)
	tIvRt3 = 1.0/math.sqrt(3.0)
	##~~~~~~~##
	tCopy = tEmpty2.copy()
	##~~~~~~~##
	tSetAllTo = 10
	##~~~~~~~##
	tResult1 = tFourVectors['a'] * tFourVectors['b']
	tAPlusB = tFourVectors['a'] + tFourVectors['b']
	tResult2 = tFourVectors['b'] * tFourVectors['a']
	tResult3 = tAPlusB * tFourVectors['c']
	##~~~~~~~##
	tTimesBy = float(random.randrange(2,5))
	tResult4 = tFourVectors['a'].copy()
	tResult4 *= tTimesBy
	tResult5 = tFourVectors['a'].copy()
	tResult5 *= -tTimesBy
	tResult6 = tFourVectors['b'].copy()
	tResult6 *= tTimesBy
	##~~~~~~~##
	tResult7 = tFourVectors['a'] + tFourVectors['b']
	tResult8 = tFourVectors['a'].copy()
	tResult8 *= -1

	##Test check_is_fourVector:##
	print "\n--------------------------------------------------\n"
	print "Testing check_is_fourVector:\n"
	tSuccessful = True
	print "Given:" , tFourVectors['a']
	print "The function returns:" , check_is_fourVector(tFourVectors['a'])
	if not check_is_fourVector(tFourVectors['a']):
		tSuccessful = False
	print "Given : " , tFloat
	print "The function returns:" , check_is_fourVector(tFloat)
	if check_is_fourVector(tFloat):
		tSuccessful = False
	print "Given : " , tString
	print "The function returns:" , check_is_fourVector(tString)
	if check_is_fourVector(tString):
		tSuccessful = False
	if tSuccessful:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_is_fourVector."
	assertions.pause(__name__)

	##Test check_different_direction:##
	print "\n--------------------------------------------------\n"
	print "Testing check_different_direction:\n"
	print "Given:" , tFourVectors['a'], "and:\n", tFourVectors['b']
	print "The function returns:" , check_different_direction(tFourVectors['a'],tFourVectors['b'])
	print "\nGiven:" , tFourVectors['a'], "and:\n", tSameDirection
	print "The function returns:" , check_different_direction(tFourVectors['a'],tSameDirection)
	print "\nGiven:" , tFourVectors['a'], "and:\n", tOppositeDirection
	print "The function returns:" , check_different_direction(tFourVectors['a'],tOppositeDirection)
	print "\nGiven:" , tFourVectors['a'], "and:\n", tOneZero
	print "The function returns:" , check_different_direction(tFourVectors['a'],tOneZero)
	print "\nGiven:" , tFourVectors['a'], "and:\n", tTwoZero
	print "The function returns:" , check_different_direction(tFourVectors['a'],tTwoZero)
	print "\nGiven:" , tFourVectors['a'], "and:\n", tAllZero
	print "The function returns:" , check_different_direction(tFourVectors['a'],tAllZero)
	print "\nExpect: True,False,True,True?,True?,False"
	print "\nFinished testing check_different_direction."
	assertions.pause(__name__)

	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "/////////////////////////"
	print "Testing fourVector class:"
	print "/////////////////////////"
	assertions.pause(__name__)

	##Test count_initially_empty, update_empty and __nonzero__:##
	print "\n--------------------------------------------------\n"
	print "Testing count_initially_empty, update_empty and __nonzero__:\n"
	print "Initallising with no variables gives:" , tEmpty1._fourVector__numberInitiallyEmpty , "empty," , tEmpty1._fourVector__positionsEmpty
	print "Initallising with 4 variable gives:" , tEmpty2._fourVector__numberInitiallyEmpty , "empty," , tEmpty2._fourVector__positionsEmpty
	print "Checking __nonzero__ returns:" , tEmpty1.__nonzero__() , "and" , tEmpty2.__nonzero__()
	if not ((tEmpty1.__nonzero__() == True) or (tEmpty2.__nonzero__() == False)):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing count_initially_empty, update_empty and __nonzero__."
	assertions.pause(__name__)

	##Test force_float_vector:##
	print "\n--------------------------------------------------\n"
	print "Testing force_float_vector:\n"
	print "Using the second vector above:" , tEmpty2
	tEmpty2.force_float_vector()
	print "gives:" , tEmpty2
	if not ((type(tEmpty2[0]) != float) or (type(tEmpty2[1]) != float) or (type(tEmpty2[2]) != float) or (type(tEmpty2[3]) != float)):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing force_float_vector."
	assertions.pause(__name__)

	##Test copy:##
	print "\n--------------------------------------------------\n"
	print "Testing copy:\n"
	print "Using the vector above gives:" , tCopy
	tCopy *= 2
	print "Multiplying the copy by 2 gives:" , tCopy
	print "whilst the original vector is:" , tEmpty2
	if tEmpty2 == tCopy / 2.0:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing copy."
	assertions.pause(__name__)

	##Test calculate_cartesian_magnitude:##
	print "\n--------------------------------------------------\n"
	print "Testing calculate_cartesian_magnitude:\n"
	print "Using:" , tMagnitude
	print "gives:" , tMagnitude.calculate_cartesian_magnitude()
	if (tMagnitude.calculate_cartesian_magnitude() == math.sqrt(5.25)):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing calculate_cartesian_magnitude."
	assertions.pause(__name__)

	##Test get_cartesian_unit_vector:##
	print "\n--------------------------------------------------\n"
	tSuccessful = True
	print "Testing get_cartesian_unit_vector:\n"
	print "Using:" , tUnitVector1
	print "which has magnitude:", tUnitVector1.calculate_cartesian_magnitude()
	print "gives:" , tUnitVector1.get_cartesian_unit_vector()
	if not (tUnitVector1.get_cartesian_unit_vector() == fourVector(0,tIvRt3,tIvRt3,tIvRt3)):
		tSuccessful = False
	print "which has magnitude:", tUnitVector1.get_cartesian_unit_vector().calculate_cartesian_magnitude()
	if not (tUnitVector1.get_cartesian_unit_vector().calculate_cartesian_magnitude() == 1.0):
		tSuccessful = False
	print "Using:" , tUnitVector2
	print "which has magnitude:", tUnitVector2.calculate_cartesian_magnitude()
	print "gives:" , tUnitVector2.get_cartesian_unit_vector()
	if not (tUnitVector2.get_cartesian_unit_vector() == fourVector(0,tIvRt3,-tIvRt3,tIvRt3)):
		tSuccessful = False
	print "which has magnitude:", tUnitVector2.get_cartesian_unit_vector().calculate_cartesian_magnitude()
	if not (tUnitVector2.get_cartesian_unit_vector().calculate_cartesian_magnitude() == 1.0):
		tSuccessful = False
	if tSuccessful:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing get_cartesian_unit_vector."
	assertions.pause(__name__)

	##Test __setitem__ and __getitem__:##
	print "\n--------------------------------------------------\n"
	print "Testing __setitem__ and __getitem__:\n"
	print "d =" , tFourVectors['d']
	for __i12 in range(4):
		print "Element" , __i12 ,"was =" , tFourVectors['d'][__i12]
		tFourVectors['d'][__i12] = tSetAllTo
	print "Set each individual element to:" , tSetAllTo
	print "d = " , tFourVectors['d']
	if ((tFourVectors['d'][0] == tSetAllTo) and (tFourVectors['d'][1] == tSetAllTo)):
		if ((tFourVectors['d'][2] == tSetAllTo) and (tFourVectors['d'][3] == tSetAllTo)):
			print "\nTest successful!"
		else:
			print "\nTest failed!"
	else:
		print "\nTest failed!"
	print "\nFinished testing __setitem__ and __getitem__."
	assertions.pause(__name__)

	##Test __str__ and __repr__:##
	print "\n--------------------------------------------------\n"
	print "Testing __str__ and __repr__:\n"
	print "Using the two test vectors from before:"
	print "__str__ gives:" , str(tEmpty1)
	print "and: " , str(tEmpty2)
	print "__repr__ gives:" , repr(tEmpty1)
	print "and: " , repr(tEmpty2)
	print "\nFinished testing __str__ and __repr__."
	assertions.pause(__name__)
	
	##Test __mul__, the Minkowski inner product:##
	print "\n---------------------------------------------\n"
	print "Testing __mul__, the Minkowski inner product:\n"
	print "a = " , tFourVectors['a']
	print "b = " , tFourVectors['b']
	print "a * b = " , tResult1
	print "b * a = " , tResult2
	print "a + b = " , tAPlusB
	print "    c = " , tFourVectors['c']
	print "(a + b) * c = " , tResult3
	expected1 = (tFourVectors['b'][0]*tFourVectors['a'][0]) - (tFourVectors['b'][1]*tFourVectors['a'][1])
	expected1 -= ((tFourVectors['b'][2]*tFourVectors['a'][2]) + (tFourVectors['b'][3]*tFourVectors['a'][3]))
	if tResult1 == expected1:
		print "\nTest (a*b) successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing four-vector __mul__, the Minkowski inner product."
	assertions.pause(__name__)

	##Test __imul__, __div__ and __idiv__:##
	tSuccessful = True
	print "\n---------------------------------------------\n"
	print "Testing __imul__, __div__ and __idiv__:\n"
	print "a = " , tFourVectors['a']
	print "a *" , tTimesBy , "=" , tResult4
	tResult4 /= tTimesBy
	print "and /" , tTimesBy , "=" , tResult4
	if not tResult4 == tFourVectors['a']:
		tSuccessful = False
	print "a *" , -tTimesBy , "=" , tResult5
	tResult5 /= -tTimesBy
	print "and /" , -tTimesBy , "=" , tResult5
	if not tResult5 == tFourVectors['a']:
		tSuccessful = False
	print "b = " , tFourVectors['b']
	print "b *" , tTimesBy , "=" , tResult6
	tResult6 /= tTimesBy
	if not tResult6 == tFourVectors['b']:
		tSuccessful = False
	print "and /" , tTimesBy , "=" , tResult6
	if tSuccessful:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing __imul__, __div__ and __idiv__."
	assertions.pause(__name__)

	##Test __add__, __iadd__, __sub__ and __isub__:##
	print "\n---------------------------------------------\n"
	print "Testing __add__, __iadd__, __sub__ and __isub__:\n"
	print "a = " , tFourVectors['a']
	print "b = " , tFourVectors['b']
	print "a + b = " , tResult7
	print "    c = " , tFourVectors['c']
	tResult7 += tFourVectors['c']
	print "(a + b) + c = " , tResult7
	tResult7 = tResult7 - tFourVectors['c']
	print "((a + b) + c) - c = " , tResult7
	tResult7 -= tFourVectors['b']
	print "(((a + b) + c) - c) - b = " , tResult7
	if tResult7 == tFourVectors['a']:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing __add__, __iadd__, __sub__ and __isub__."
	assertions.pause(__name__)

	##Test __abs__:##
	print "\n---------------------------------------------\n"
	tSuccessful = True
	print "Testing __abs__:\n"
	print "a = " , tFourVectors['a']
	print "and has absolute value = ", abs(tFourVectors['a'])
	print "b = " , tFourVectors['b']
	print "and has absolute value = ", abs(tFourVectors['b'])
	print "-a = " , tResult8
	print "and has absolute value = ", abs(tResult8)
	if not abs(tFourVectors['a']) == tFourVectors['a']*tFourVectors['a']:
		tSuccessful = False
	if not abs(tResult8) == abs(tFourVectors['a']):
		tSuccessful = False
	if tSuccessful:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing __abs__."
	assertions.pause(__name__)

	##Test __eq__ and __ne__, the compairson operators:##
	print "\n---------------------------------------------\n"
	print "Testing __eq__ and __ne__, the comparison operators:\n"
	tSuccessful = True
	print "a = " , tFourVectors['a']
	print "b = " , tFourVectors['b']
	print "a == b?" , tFourVectors['a'] ==  tFourVectors['b']
	if tFourVectors['a'] ==  tFourVectors['b']:
		tSuccessful = False
	print "a != b?" , tFourVectors['a'] != tFourVectors['b']
	if not tFourVectors['a'] !=  tFourVectors['b']:
		tSuccessful = False
	print "a == a?" , tFourVectors['a'] ==  tFourVectors['a']
	if not tFourVectors['a'] ==  tFourVectors['a']:
		tSuccessful = False
	print "a != a?" , tFourVectors['a'] != tFourVectors['a']
	if tFourVectors['a'] !=  tFourVectors['a']:
		tSuccessful = False
	print "b == b?" , tFourVectors['b'] ==  tFourVectors['b']
	if not tFourVectors['b'] ==  tFourVectors['b']:
		tSuccessful = False
	print "b != b?" , tFourVectors['b'] != tFourVectors['b']
	if tFourVectors['b'] !=  tFourVectors['b']:
		tSuccessful = False
	if tSuccessful:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing __eq__ and __ne__, the comparison operators."

	##Done testing:##
	print "\n---------------------------------------------\n"
	print "/////////////////////////////////////"
	print "Finished checking fourVectors module!"
	print "/////////////////////////////////////"