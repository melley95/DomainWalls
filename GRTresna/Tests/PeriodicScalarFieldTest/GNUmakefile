# -*- Mode: Makefile -*- 

# the location of the Chombo "lib" directory
ifndef CHOMBO_HOME
    $(error Please define CHOMBO_HOME - see installation instructions.)
endif

# trace the chain of included makefiles
makefiles += releasedExamples_AMRPoisson_execVariableCoefficient

# the base name(s) of the application(s) in this directory
ebase = PeriodicScalarTest

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRElliptic AMRTools BoxTools

# input file for 'run' target
INPUT = params.txt

# application-specific targets
src_dirs := ../../Source \
            ../../Source/Core \
            ../../Source/Matter \
            ../../Source/Methods \
            ../../Source/Tools \
            ../../Source/Variables \
            ../../Source/TaggingCriteria \
  	        ../../Source/Operator \
            ../../Source/Operator/SolverOperator 

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.test
