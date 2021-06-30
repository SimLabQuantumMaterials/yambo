define mk_yambo
 LLIBS="";for exe in $$XLIBS; do LLIBS="$$LLIBS -l$$exe" ; done ; \
 for exe in $$X2DO; do \
  DLIBS="-ldriver";for lib in $$DRILIBS; do DLIBS="$$DLIBS -l$$exe$$lib" ; done ; \
  if test ! -f "$(bindir)/$$exe" || test "$(keep_objs)" = yes  ; then \
   if test ! -d $$XPATH ; then mkdir -p $$XPATH ; fi ; \
   if [ "$(compdir)" != "$(prefix)" ] && [ -f $$VPATH/.objects ] ; then \
    cp $$VPATH/.objects $$XPATH ; \
   fi ; \
   ./sbin/compilation/helper.sh -d $$XPATH -t $$exe -o .objects -m x -- "$$DLIBS $$LLIBS $(xcpp) $$ADF"; \
   cd $$XPATH ; $(MAKE) VPATH=$$VPATH exe || exit "$$?" ; \
  fi ; \
  echo " " ; \
 done
endef
define mk_ypp
 LLIBS="";for exe in $$XLIBS; do LLIBS="$$LLIBS -l$$exe" ; done ; \
 for exe in $$X_ypp_LIBS; do LLIBS="$$LLIBS -l_ypp_$$exe" ; done ; \
 for exe in $$X2DO; do \
  DLIBS="-ldriver";for lib in $$DRILIBS; do DLIBS="$$DLIBS -l$$exe$$lib" ; done ; \
  if test ! -f "$(bindir)/$$exe" || test "$(keep_objs)" = yes  ; then \
   if [ "$(compdir)" != "$(prefix)" ] && [ -f $$VPATH/.objects ] ; then \
    cp $$VPATH/.objects $$XPATH ; \
   fi ; \
   ./sbin/compilation/helper.sh -d $$XPATH -t $$exe -o .objects -m x -- "$$DLIBS $$LLIBS $(xcpp) $$ADF" ; \
   cd $$XPATH ; $(MAKE) VPATH=$$VPATH exe || exit "$$?" ; \
  fi ; \
  echo " " ; \
 done
endef
