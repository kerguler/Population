AM_CFLAGS = -I$(top_srcdir)/src/lib -Wall

AM_CFLAGS += $(GLIB_CFLAGS)

# The libraries to build
lib_LTLIBRARIES =  libpopulation.la

# The header files to include
libpopulation_ladir = $(includedir)/population

# The files to add to the library and to the source distribution
# ---
libpopulation_la_HEADERS = uthash.h population.h

libpopulation_la_SOURCES = \
	                 $(libpopulation_la_HEADERS) \
				     ran_gen.c \
				     gamma.c \
					 population.c

libpopulation_la_LIBADD = $(GLIB_LIBS)
