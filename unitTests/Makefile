# Makefile for unit tests

test: test.all

# ----- All unit tests go below ------------------

unit.NameType:
	@-cd NameType && ./UnitTest.sh $(OPT)

unit.ArgList:
	@-cd ArgList && ./UnitTest.sh $(OPT)

# ----- Every unit test should go here -----------
COMPLETETESTS= \
  unit.NameType \
  unit.ArgList

test.cpptraj: $(COMPLETETESTS)

# ------------------------------------------------

test.complete: CpptrajTest.sh UnitMaster.sh
	@./CpptrajTest.sh --target test.cpptraj $(OPT)

test.all:
	$(MAKE) test.complete summary

test.showerrors:
	$(MAKE) test.complete summary OPT="$(OPT) showerrors"

summary: CpptrajTest.sh UnitMaster.sh 
	@./CpptrajTest.sh summary $(OPT)

clean:
	@echo "Cleaning unit tests."
	$(MAKE) test.complete OPT=clean
