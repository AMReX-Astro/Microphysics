# This script takes a list of ReacLib version 1 format rate files and creates
# a Fortran 90 module, rates_module, which contains storage for common
# temperature factors, a routine to initialize those factors, and then
# a routine that calculates the reaction rate for each reaction in the rate
# files.
#
# This is useful for quickly building a nuclear reaction network.  One still
# must specify how these rates are to be used in the RHS of the ODEs.
#
# CMM: March 2012
# CMM: June  2012 --- added temperature derivative support
#
import sys

# formatting information about ReacLib v. 1 format
numHeaderLines = 3
lengthEntry = 13
numCharLabels = 35
numTerms = 7
zeroString = "0.0"

def modulePreamble():
    """Start the rate_module by declaring storage for common temperature
    factors and any possible weak rates.  Return this text"""
    text="""
! This module is used to calculate reaction rates.
! Storage is used to handle common temperature factors used
! in the various reaction rates.
!
! This module is created via the make_rate_module.py routine.
module rates_module

  use bl_types
  use bl_constants_module

  implicit none

  type temp_t
     real (kind=dp_t) :: t9
     real (kind=dp_t) :: t9i
     real (kind=dp_t) :: t9i13
     real (kind=dp_t) :: t913
     real (kind=dp_t) :: t953
     real (kind=dp_t) :: lnt9
  end type temp_t

"""
    return text

def tFactorData():
    """Add  a routine to initialize the tfactors given the temperature.
    Return this text."""
    
    text="""

contains

  function calc_tfactors(t9) result (tfactors)

    real (kind=dp_t), intent(in   ) :: t9
    type (temp_t) :: tfactors

    tfactors%t9 = t9
    tfactors%t9i = 1.d0 / tfactors%t9
    tfactors%t9i13 = tfactors%t9i**THIRD
    tfactors%t913 = tfactors%t9**THIRD
    tfactors%t953 = tfactors%t9 * tfactors%t913 * tfactors%t913
    tfactors%lnt9 = log(tfactors%t9)

  end function calc_tfactors

"""
    return text

def buildRateRoutine(name,aFactors):
    """Given a reaction name and the list of aFactors from the ReacLib rate
       file, build the appropriate structure for the reaction rate as a 
       subroutine.  Return this text."""
    text= """
  subroutine %s(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    real (kind=dp_t), intent(  out) :: rate
    real (kind=dp_t), intent(  out) :: dratedt

    real (kind=dp_t) :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

""" % name

    nsets = len(aFactors)/numTerms
    vString =[]
    dvString =[]
    facString = []
    rString = ["    rate =",]
    drString = ["    dratedt =",]
    for set in range(nsets):
        start = set*numTerms
        vString.append("r%s" % set)
        dvString.append("dr%sdt" % set)

        str = ""
        if not aFactors[start+1] == 0.0:
            str += "    ct9i = tfactors%%t9i*(%sd0)\n" % (aFactors[start+1])
        if not aFactors[start+2] == 0.0:
            str += "    ct9i13 = tfactors%%t9i13*(%sd0)\n" % (aFactors[start+2])
        if not aFactors[start+3] == 0.0:
            str += "    ct913 = tfactors%%t913*(%sd0)\n" % (aFactors[start+3])
        if not aFactors[start+4] == 0.0:
            str += "    ct9 = tfactors%%t9*(%sd0)\n" % (aFactors[start+4])
        if not aFactors[start+5] == 0.0:
            str += "    ct953 = tfactors%%t953*(%sd0)\n" % (aFactors[start+5])
        if not aFactors[start+6] == 0.0:
            str += "    clnt9 = tfactors%%lnt9*(%sd0)\n" % (aFactors[start+6])

        facString.append((str))

        str = "    %s = exp(%sd0 &\n" % (vString[set], aFactors[start])
        if not aFactors[start+1] == 0.0:
            str += "         +ct9i &\n"
        if not aFactors[start+2] == 0.0:
            str += "         +ct9i13 &\n"
        if not aFactors[start+3] == 0.0:
            str += "         +ct913 &\n"
        if not aFactors[start+4] == 0.0:
            str += "         +ct9 &\n"
        if not aFactors[start+5] == 0.0:
            str += "         +ct953 &\n"
        if not aFactors[start+6] == 0.0:
            str += "         +clnt9)"
        else:
            str += "         )"

        facString.append((str))

        str = "    %s = %s * tfactors%%t9i * ( &\n" % (dvString[set],vString[set])
        if not aFactors[start+1] == 0.0:
            str += "         -ct9i &\n"
        if not aFactors[start+2] == 0.0:
            str += "         -THIRD*ct9i13 &\n"
        if not aFactors[start+3] == 0.0:
            str += "         +THIRD*ct913 &\n"
        if not aFactors[start+4] == 0.0:
            str += "         +ct9 &\n"
        if not aFactors[start+5] == 0.0:
            str += "         +FIVE3RD*ct953 &\n" 
        if not aFactors[start+6] == 0.0:
            str += "         +(%sd0))" % (aFactors[start+6])
        else:
            str += "         )"

        facString.append((str))
        facString.append("\n")

        rString.append("%s +" % vString[set])
        drString.append("%s +" % dvString[set])

    vString = ", ".join(vString)
    dvString = ", ".join(dvString)
    facString = "\n".join(facString)
    rString = " ".join(rString); rString = rString[:-1]
    drString = " ".join(drString); drString = drString[:-1]

    text += ("    real (kind=dp_t) :: %s\n"
             +"    real (kind=dp_t) :: %s\n") % (vString,dvString)
    text +="""
    rate = ZERO
    dratedt = ZERO
"""
    text += "\n%s\n" % facString
    text += "\n%s\n%s" % (rString,drString)

    text += """

  end subroutine %s

""" % name
    return text

def numReactantsProducts(chapter):
    """Return the number of reactants and products (in particles) for a 
    reaction of ReacLib Chapter type chapter."""
    r,p = 1,1
    if chapter in [2,5,9,10]: p+= 1
    if chapter in [3,6]: p+= 2
    if chapter in [7,11]: p+= 3
    if chapter > 3 and chapter is not 11: r+= 1
    if chapter > 7 and chapter is not 11: r+= 1
    if chapter is 10: r+= 1

    return r,p

def parseRateFile(file):
    """Given a ReacLib v 1 rate file, file, parse it and determine an 
    appropriate name and the "a" factors that build the rate.  Return the
    name as a string and the "a" factors as a list."""

    try: 
        fh = open(file,'r')
    except IOError:
        print "Couldn't open file %f for reading" % file
        sys.exit()

    data = fh.readlines()
    fh.close()

    r,p = numReactantsProducts(int(data[0]))

    foundLabels = False
    rateString = "rate"
    
    aFactors = []

    for line in data[numHeaderLines:]:
        if line.startswith('  '): # this is a header line

            if not foundLabels:   # only do this once
                labels = line[:numCharLabels].split()

                for i,label in enumerate(labels):
                    if i == r: rateString += "_to"
                    rateString += "_%s" % label
                foundLabels = True
        else:
            # this is a data line; some have 4 entries, some have 3
            nterms = 4 if len(line.strip()) > 3*lengthEntry else 3
            for term in range(nterms):
                start = term*lengthEntry
                stop = start + lengthEntry
                
                aFactors.append(float(line[start:stop]))

    if len(aFactors) % numTerms is not 0:
        print "Error parsing rate file %s: missing aFactors!" % file
        sys.exit()

    return rateString,aFactors

def rateIsWeakDecay(factors):
    return all([repr(f) == zeroString for f in factors[1:]])

def weakRate(name,factor):

    textParameter="""  real (kind=dp_t), parameter :: %s = exp(%sd0)
""" % (name[5:], factor)

    textSubroutine="""
  subroutine %s(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    real (kind=dp_t), intent(  out) :: rate
    real (kind=dp_t), intent(  out) :: dratedt


    rate = %s
    dratedt = ZERO

  end subroutine

""" % (name,name[5:])

    return textParameter, textSubroutine

def make_rates_module(rateFiles):
    """Build the F90 module containing the subroutines for the reactions 
    given in the list of ReacLib rate files, rateFiles, along with the helper 
    routine to set the common temperature factors.  Return all this as a string
    for easy piping."""

    # the rates_module contains the rate routines and the tfactors
    preamble = modulePreamble()

    # build the tfactors data routine; this is the start of
    # functions/subroutines
    out_text = tFactorData()

    # add the rate routines
    for file in rateFiles:
        name, factors = parseRateFile(file)
        if (rateIsWeakDecay(factors)):
            parameter, subroutine = weakRate(name,factors[0])
            preamble += parameter
            out_text += subroutine
        else:
            out_text += buildRateRoutine(name,factors)

    # close rates_module
    out_text += """

end module rates_module"""

    return preamble + out_text


if __name__ == "__main__":
    # dump the output module to screen
    print make_rates_module(sys.argv[1:])
