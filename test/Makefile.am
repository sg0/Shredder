check_PROGRAMS = app create read stars cliques

app_CFLAGS = -I$(top_srcdir)/src
create_CFLAGS = -I$(top_srcdir)/src
read_CFLAGS = -I$(top_srcdir)/src
stars_CFLAGS = -I$(top_srcdir)/src
cliques_CFLAGS = -I$(top_srcdir)/src

app_SOURCES = example.c
create_SOURCES = create.c
read_SOURCES = read.c
stars_SOURCES = stars.c
cliques_SOURCES = cliques.c

nodist_app_SOURCES = $(app_src_SOURCES)
nodist_create_SOURCES = $(app_src_SOURCES)
nodist_read_SOURCES = $(app_src_SOURCES)
nodist_stars_SOURCES = $(app_src_SOURCES)
nodist_cliques_SOURCES = $(app_src_SOURCES)

app_src_SOURCES = shredder.c

CLEANFILES = $(app_src_SOURCES)

$(app_src_SOURCES):
	echo "#include \"$(top_builddir)/src/$@\"" >$@