## w32bc5.mak (W32 single-threaded static library, BC++ 5.2) ##

##----------------------------------------------------------------------
## To build W32 single-threaded static library for GLPK with BC++ 5.2
## enter the subdirectory where the file you are reading is placed and
## type the command:
##
##    make.exe -f w32bc5.mak
##
## The following two files will be created in the same subdirectory:
##
##    glpk.lib,   GLPK static library, and
##
##    glpsol.exe, stand-alone LP/MIP solver.
##
## To check if GLPK has been successfully compiled type the command:
##
##    make.exe -f w32bc5.mak check
##
## Header files needed to use the GLPK static library are contained in
## the subdirectory 'include'.
##
## Written by Andrew Makhorin <mao@mai2.rcnet.ru>, <mao@gnu.org>.
##----------------------------------------------------------------------

CFLAGS = -wall

.c.obj:
        bcc32.exe $(CFLAGS) -Iinclude -o$*.obj -c $*.c

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
        sysdep\w32\glplib1b.obj \
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
        tlib.exe @&&|
        glpk.lib /C /0 +$(**:.obj=.obj+)
|

glpsol.exe: examples\glpsol.obj glpk.lib
        bcc32.exe $(CFLAGS) -eglpsol.exe examples\glpsol.obj glpk.lib

tspsol.exe: examples\tspsol.obj glpk.lib
        bcc32.exe $(CFLAGS) -etspsol.exe examples\tspsol.obj glpk.lib

## eof ##
