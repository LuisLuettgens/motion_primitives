############################## 
# Makefile-Vorlage für Latex #     
##############################

MAINTEXFILE := TransWORHP.tex

##################################################
# Change from here ...
#
## MAINTEXFILE := test1 test2 ...
#
# ... until here
##################################################

ALLTEXFILES := $(wildcard *.tex) mk.sty 

MAINTEXFILE := $(subst .tex,,$(MAINTEXFILE))

PSFILES = $(MAINTEXFILE:=.ps)
PDFFILES = $(MAINTEXFILE:=.pdf)
DVIFILES = $(MAINTEXFILE:=.dvi)

# What output do you want?
all: $(PSFILES) 
#$(PDFFILES) 

# Auto re-running latex!
%.dvi: %.tex $(ALLTEXFILES)
	latex $<
	@logfile='$(subst tex,log,$<)' ; \
	latex_count=3 ; \
	while egrep -s 'Rerun to get cross-references right' $$logfile && [ $$latex_count -gt 0 ] ;\
	    do \
		  echo ;\
		  echo ;\
		  echo ;\
	      echo "Rerunning latex...." ;\
	      latex $< ;\
	      latex_count=`expr $$latex_count - 1` ;\
	    done


%.ps: %.dvi
	dvips $< -o $@

# Auto re-running pdflatex!
%.pdf: %.tex $(ALLTEXFILES)
	pdflatex $<
	@logfile='$(subst tex,log,$<)' ; \
	latex_count=3 ; \
	while egrep -s 'Rerun to get cross-references right' $$logfile && [ $$latex_count -gt 0 ] ;\
	    do \
		  echo ;\
		  echo ;\
		  echo ;\
	      echo "Rerunning latex...." ;\
	      pdflatex $< ;\
	      latex_count=`expr $$latex_count - 1` ;\
	    done

%: %.ps


clean:
	rm -f *.aux *.log *.tex~ *.idx *.ilg *.ind *.lol *.out *.toc *.lot *.lof
	cd programme && make clean

clean2:
	rm -f $(PSFILES) $(PDFFILES) $(DVIFILES)

%.pic: %.flw
	flow <$< >$@

bib:
#	bibtex link
	bibtex doku


matthias.bst: matthias.dbj
	latex matthias.dbj

ref:
	texref -unref skript


SPELLTEX := $(subst .tex,.spell,$(ALLTEXFILES))

spell: $(SPELLTEX)
%.spell: %.tex
	aspell --lang=de_DE-neu -c $<

ispell: skript.tex
	ispell -t -b -d german -ktexskip2 +psfrag,fref,rnode,pscircle,begin,psaxes,psframe,psline,ncline,cnodeput,citeyear -ktexskip1 +fref,psset,citet,citeauthor $<

