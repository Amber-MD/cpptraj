# Amber buildfile for cpptraj's internal readline.
# This readline has been modified to not depend on an external termcap library.

set(CPPTRAJ_READLINE_SOURCES
	bind.c 
	callback.c 
	compat.c 
	complete.c 
	display.c 
	funmap.c 
	histexpand.c 
	histfile.c 
	history.c 
	histsearch.c 
	input.c 
	isearch.c 
	keymaps.c 
	kill.c 
	macro.c 
	mbutil.c 
	misc.c 
	nls.c 
	parens.c 
	readline.c 
	rltty.c 
	savestring.c 
	search.c 
	shell.c 
	signals.c 
	termcap.c 
	terminal.c 
	text.c 
	tilde.c 
	tparam.c 
	undo.c 
	util.c 
	version.c 
	vi_mode.c 
	xmalloc.c)

add_compile_options(${OPT_CFLAGS})

add_library(readline STATIC ${CPPTRAJ_READLINE_SOURCES})
target_compile_definitions(readline PRIVATE HAVE_CONFIG_H=1) # make sure the code uses the premade config.h
make_pic_if_needed(readline)
target_include_directories(readline PUBLIC .)
