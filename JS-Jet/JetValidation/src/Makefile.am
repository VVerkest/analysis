AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  JetValidation.h \
  jetBackgroundCut.h

ROOTDICTS = 

pcmdir = $(libdir)
nobase_dist_pcm_DATA =

libjetbackgroundcut_io_la_SOURCES = \
  $(ROOTDICTS) 

libjetbackgroundcut_la_SOURCES = \
  jetBackgroundCut.cc

libjetbackgroundcut_io_la_LIBADD = \
  -lphool -lSubsysReco -lmbd_io -ljetbase

libjetbackgroundcut_la_LIBADD = \
  libjetbackgroundcut_io.la \
  -lcalo_io \
  -lSubsysReco -lffaobjects -lmbd_io \
  -ljetbase -lglobalvertex

lib_LTLIBRARIES = \
  libJetValidation.la \
  libjetbackgroundcut_io.la \
  libjetbackgroundcut.la

libJetValidation_la_SOURCES = \
  JetValidation.cc \
  $(ROOTDICTS)

libJetValidation_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lcalotrigger \
  -lcalotrigger_io \
  -ljetbackgroundcut \
  -lfun4all \
  -lg4dst \
  -lg4eval \
  -lg4jets

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libJetValidation.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
