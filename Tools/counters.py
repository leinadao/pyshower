####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module for creating and handling counters to use during dipole showering."""

##Import required modules:##

print "\n////////////////////////"
print "Loading counters module:"
print "////////////////////////\n"

##Functions:##

def check_is_counter(toCheck):
	"""A function to check for an instance of the counter class."""
	return isinstance(toCheck,counter)

##Classes:##

class counter(object):
	"""A class for counting in dipole showering."""

	def __init__(self,startValue):
		assert (type(startValue) == int)
		self.__startValue = startValue
		self.__nextValue = self.__startValue

	def next(self):
		"""A function to return the next counter value."""
		__next = self.__nextValue
		self.__nextValue += 1
		return __next

	def count(self):
		"""A function to count the next counter value without returning it."""
		self.__nextValue += 1

	def counted(self):
		"""A function to return how many values have been counted by the counter."""
		return self.__nextValue - self.__startValue

##Set up counters:##
gluonProdCounter = counter(1)
gluonSplitCounter = counter(1)
photonProdCounter = counter(1)
kPerpProdWarningCounter = counter(1)

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##
	import assertions
	import fourVectors

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "////////////////////////"
	print "Testing counters module:"
	print "////////////////////////"
	assertions.pause(__name__)
	
	##Setup here:##
	tCounter1 = counter(0)
	tCounter2 = counter(501)
	tCounter3 = counter(0)
	tCounter4 = counter(501)
	tNumIts = 1000000

	##Test check_is_counter:##
	print "\n--------------------------------------------------\n"
	print "Testing check_is_counter:"
	print "Calling check_is_counter on instance: " , check_is_counter(tCounter1)
	print "Calling check_is_counter on wrong instance: " , check_is_counter(fourVectors.fourVector())
	print "Calling check_is_counter on 1.055: " , check_is_counter(1.055)
	print "Calling check_is_counter on 'word': " , check_is_counter('word')
	tResults = [check_is_counter(tCounter1),check_is_counter(fourVectors.fourVector())]
	tResults += [check_is_counter(1.055),check_is_counter('word')]
	if (sum(tResults) == 1):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_is_counter."
	assertions.pause(__name__)

	##Test counter class:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "//////////////////////"
	print "Testing counter class:"
	print "//////////////////////"
	assertions.pause(__name__)

	##__init__() tested implicitly.

	##Test next and counted functions:##
	print "\n--------------------------------------------------\n"
	print "Testing next and counted functions:\n"
	print "Using test counter 1 which starts at 0:"
	for ti in range(10):
		print "Calling next returns:", tCounter1.next()
	print "Calling counted returns:", tCounter1.counted()
	if ((tCounter1.counted() == 10) and (tCounter1.next() == 10)):
		print "Test successful!"
	else:
		print "Test failed!"
	assertions.pause(__name__)
	print "Using test counter 2 which starts at 501:"
	for ti in range(20):
		print "Calling next returns:", tCounter2.next()
	print "Calling counted returns:", tCounter2.counted()
	print "Call " +  str(tNumIts) + " more times:"
	tRuni = 0
	while tRuni < tNumIts:
		tA = tCounter2.next() ##Use a = to prevent print slowing test down.
		tRuni += 1
	print "Finally counted returns:", tCounter2.counted()
	if ((tCounter2.counted() == tNumIts+20) and (tCounter2.next() == tNumIts+20+501)):
		print "Test successful!"
	else:
		print "Test failed!"
	print "\nFinished testing next and counted functions."
	assertions.pause(__name__)

	##Test count and counted functions:##
	print "\n--------------------------------------------------\n"
	print "Testing count and counted functions:\n"
	print "Using test counter 3 which starts at 0:"
	for ti in range(10):
		print "Calling count:"
		tCounter3.count()
	print "Calling counted returns:", tCounter3.counted()
	if ((tCounter3.counted() == 10) and (tCounter3.next() == 10)):
		print "Test successful!"
	else:
		print "Test failed!"
	assertions.pause(__name__)
	print "Using test counter 4 which starts at 501:"
	for ti in range(20):
		print "Calling count:"
		tCounter4.count()
	print "Calling counted returns:", tCounter4.counted()
	print "Call " +  str(tNumIts) + " more times:"
	tRuni = 0
	while tRuni < tNumIts:
		tCounter4.count()
		tRuni += 1
	print "Finally counted returns:", tCounter4.counted()
	if ((tCounter4.counted() == tNumIts+20) and (tCounter4.next() == tNumIts+20+501)):
		print "Test successful!"
	else:
		print "Test failed!"
	print "\nFinished testing count and counted functions."
	assertions.pause(__name__)

	##Done testing:##
	print "\n---------------------------------------------\n"
	print "//////////////////////////////////"
	print "Finished checking counters module!"
	print "//////////////////////////////////"