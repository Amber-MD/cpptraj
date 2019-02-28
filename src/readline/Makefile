# Makefile for readline bundled with cpptraj.
include readline.config.h

AR = ar cqs

SOURCES       = bind.c \
		callback.c \
		compat.c \
		complete.c \
		display.c \
		funmap.c \
		histexpand.c \
		histfile.c \
		history.c \
		histsearch.c \
		input.c \
		isearch.c \
		keymaps.c \
		kill.c \
		macro.c \
		mbutil.c \
		misc.c \
		nls.c \
		parens.c \
		readline.c \
		rltty.c \
		savestring.c \
		search.c \
		shell.c \
		signals.c \
		termcap.c \
		terminal.c \
		text.c \
		tilde.c \
		tparam.c \
		undo.c \
		util.c \
		version.c \
		vi_mode.c \
		xmalloc.c

OBJECTS=$(SOURCES:.c=.o)

all: libreadline.a

libreadline.a: $(OBJECTS)
	$(AR) libreadline.a $(OBJECTS)

uninstall: clean

clean:
	/bin/rm -f *.o libreadline.a

# Dependency targets
../findDepend: ../FindDepend.o 
	$(CXX) -o ../findDepend ../FindDepend.o

depend: ../findDepend
	../findDepend $(SOURCES) > readlinedepend

dependclean:
	/bin/rm -f ../FindDepend.o ../findDepend

# Dependencies
include readlinedepend
