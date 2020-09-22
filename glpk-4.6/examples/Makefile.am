## examples/Makefile.am ##

INCLUDES = -I../include

LDADD = ../src/libglpk.a

bin_PROGRAMS = glpsol tspsol

glpsol_SOURCES = glpsol.c

tspsol_SOURCES = tspsol.c

check: glpsol
	./glpsol ./plan.mps
	@echo ""
	@echo "If optimal solution has been found, the check is OK"
	@echo ""

## eof ##
