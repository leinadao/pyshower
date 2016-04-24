####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module for creating data loggers to use during dipole showering."""

##Import required modules:##

print "\n///////////////////////////"
print "Loading dataLoggers module:"
print "///////////////////////////\n"

##Functions:##

def check_is_dataLogger(toCheck):
	"""A function to check for an instance of the dataLogger class."""
	return isinstance(toCheck,dataLogger)

##Classes:##

class dataLogger(object):
	"""A class for storing results from dipole showering."""
	def __init__(self):
		self.__logged = []

	def store(self,toLog):
		"""A function to log the next value."""
		self.__logged.append(toLog)

	def output(self):
		"""A function to return how many values have been issued by the dataLogger."""
		return self.__logged

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##
	import random
	import assertions
	import fourVectors

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "///////////////////////////"
	print "Testing dataLoggers module:"
	print "///////////////////////////"
	assertions.pause(__name__)
	
	##Setup here:##
	tLogger1 = dataLogger()
	tNumIts = 1000000

	##Test check_is_dataLogger:##
	print "\n--------------------------------------------------\n"
	print "Testing check_is_dataLogger:"
	print "Calling check_is_dataLogger on instance: " , check_is_dataLogger(tLogger1)
	print "Calling check_is_dataLogger on wrong instance: " , check_is_dataLogger(fourVectors.fourVector())
	print "Calling check_is_dataLogger on 1.055: " , check_is_dataLogger(1.055)
	print "Calling check_is_dataLogger on 'word': " , check_is_dataLogger('word')
	tResults = [check_is_dataLogger(tLogger1),check_is_dataLogger(fourVectors.fourVector())]
	tResults += [check_is_dataLogger(1.055),check_is_dataLogger('word')]
	if (sum(tResults) == 1):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_is_dataLogger."
	assertions.pause(__name__)

	##Test dataLogger class:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "/////////////////////////"
	print "Testing dataLogger class:"
	print "/////////////////////////"
	assertions.pause(__name__)

	##__init__() tested implicitly.

	##Test store and output functions:##
	print "\n--------------------------------------------------\n"
	print "Testing store and output functions:\n"
	print "Using test dataLogger1 which starts at 0:"
	print "Storing", tNumIts, "random numbers:"
	tL = []
	for ti in range(tNumIts):
		tR = random.random()
		tLogger1.store(tR)
		tL.append(tR)
	print "Testing output list matches that expected:\n"
	if (tLogger1.output() == tL):
		print "Test successful!"
	else:
		print "Test failed!"
	print "\nFinished testing store and output functions."
	assertions.pause(__name__)

	##Done testing:##
	print "\n---------------------------------------------\n"
	print "/////////////////////////////////////"
	print "Finished checking dataLoggers module!"
	print "/////////////////////////////////////"