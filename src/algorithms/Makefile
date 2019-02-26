##                           -*- Mode: Makefile -*- 
## Makefile --- 
## Author          : Truong Nguyen and Marcel Turcotte
## Created On      : Wed Jun  8 10:53:00 2005
## Last Modified By: Marcel Turcotte
## Last Modified On: Mon Aug  8 18:06:00 2005
##
## This copyrighted source code is freely distributed under the terms
## of the GNU General Public License. 
## See the files COPYRIGHT and LICENSE for details.

SHELL = /bin/sh

OBJECTS = seed.o ida.o stems.o motif.o misc.o
SEED_W_OBJECTS = seed_w.o match_w.o kfunc.o

BINARIES = seed find match
#  match_w

LIBS = -lbio -lvtree -ldev -lm
LIBDIR = -L../libbio -L../libvtree -L../libdev
INCDIR = -I../libbio -I../libvtree -I../libdev

.SUFFIXES:
.SUFFIXES: .c .o

.c.o: 
	$(CC) $(CFLAGS) $(RNALIB) $(INCDIR) $(RNALIB_INCLUDE) -c $<

all: $(BINARIES)

seed: $(OBJECTS)
	$(CC) -o seed $(OBJECTS) $(CFLAGS) $(LIBDIR) $(RNALIB_LIB) $(LDFLAGS) $(LIBS) $(RNALIB_LIBS)

seed_w: $(SEED_W_OBJECTS)
	$(CC) -o seed_w $(SEED_W_OBJECTS) $(CFLAGS) $(LIBDIR) $(LIBS)

kfunc: kfunc.o
	$(CC) -o kfunc kfunc.o $(LIBS)

match_w: match_w.o
	$(CC) -o match_w match_w.o $(CFLAGS) $(LIBDIR) $(LDFLAGS) $(LIBS)

find:  find.o
	$(CC) -o find find.o $(CFLAGS) $(LIBDIR) $(LDFLAGS) $(LIBS)

match:  match.o
	$(CC) -o match match.o $(CFLAGS) $(LIBDIR) $(LDFLAGS) $(LIBS)

install: $(BINARIES)
	cp $(BINARIES) $(BIN_DIR)

check:
	true

clean:
	rm -f $(OBJECTS) $(SEED_W_OBJECTS) *~ seed seed_w test find find.o match match.o match_w match_w.o kfunc kfunc.o
