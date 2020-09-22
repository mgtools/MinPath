## d32dmc.mak ##

##----------------------------------------------------------------------
## Build DOS32 version of GLPK
## Needs Digital Mars C compiler from http://www.digitalmars.com
## Needs 32 bit DOS extender from http://www.dosextender.com
##
## Written by Vasily Zatsepin <zwb@online.ru>.
##
## To compile GLPK with Digital Mars C for DOS32 enter the subdirectory
## where the file you are reading is placed and type the command:
##
##    make.exe -f d32dmc.mak
##
## The following two files will be created in the same subdirectory:
##
##    glpk.lib,   GLPK static library, and
##
##    glpsol.exe, stand-alone LP/MIP solver.
##
## To check if GLPK has been successfully compiled type the command:
##
##    make.exe -f d32dmc.mak check
##
## Header files needed to use the GLPK static library are contained in
## the subdirectory 'include'.
##----------------------------------------------------------------------

CFLAGS = -mx -w

.c.obj:
        sc.exe $(CFLAGS) -Iinclude -o$*.obj -c $*.c

all: glpk.lib glpsol.exe tspsol.exe

check: glpsol.exe
      glpsol.exe examples\plan.mps

glpk.lib: src\glpavl.obj \
        src\glpchol.obj \
        src\glpdmp.obj \
        src\glpiet.obj \
        src\glpinv.obj \
        src\glpios1.obj \
        src\glpios2.obj \
        src\glpios3.obj \
        src\glpipm.obj \
        src\glplib1a.obj \
        src\glplib1b.obj \
        src\glplib2.obj \
        src\glplib3.obj \
        src\glplpp1.obj \
        src\glplpp2.obj \
        src\glplpt.obj \
        src\glplpx1.obj \
        src\glplpx2.obj \
        src\glplpx3.obj \
        src\glplpx4.obj \
        src\glplpx5.obj \
        src\glplpx6a.obj \
        src\glplpx6b.obj \
        src\glplpx6c.obj \
        src\glplpx6d.obj \
        src\glplpx7.obj \
        src\glplpx8a.obj \
        src\glplpx8b.obj \
        src\glplpx8c.obj \
        src\glplpx8d.obj \
        src\glplpx8e.obj \
        src\glpluf.obj \
        src\glpmat.obj \
        src\glpmip1.obj \
        src\glpmip2.obj \
        src\glpmpl1.obj \
        src\glpmpl2.obj \
        src\glpmpl3.obj \
        src\glpmpl4.obj \
        src\glpmps.obj \
        src\glpqmd.obj \
        src\glprng.obj \
        src\glpspx1.obj \
        src\glpspx2.obj \
        src\glpstr.obj \
        src\glptsp.obj
        lib.exe -c $@ $**

glpsol.exe: examples\glpsol.obj glpk.lib
        sc.exe $(CFLAGS) examples\glpsol.obj glpk.lib X32.lib

tspsol.exe: examples\tspsol.obj glpk.lib
        sc.exe $(CFLAGS) examples\tspsol.obj glpk.lib X32.lib

## eof ##
