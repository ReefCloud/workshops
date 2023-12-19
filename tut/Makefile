TUT_SRC_DIR  = $(addprefix , tut)
ifdef FILE
TUT_FILES = $(TUT_SRC_DIR)/$(FILE)
else
TUT_FILES ?= $(foreach sdir, $(TUT_SRC_DIR), $(wildcard $(sdir)/*.qmd))
endif

HTML_FILES := $(patsubst %.qmd, %.html, $(TUT_FILES))
PUBLIC_DIR := public/
DOCS_DIR   := $(addprefix ,docs/)

##DOCS_DIR := docs/

$(info ************************************)
$(info TUT Source directory:     $(TUT_SRC_DIR))
$(info TUT Source files:         $(TUT_FILES))
$(info TUT HTML files:           $(HTML_FILES))
$(info PUBLIC_DIR:              $(PUBLIC_DIR))
$(info DOCS_DIR:              $(DOCS_DIR))
$(info ************************************)

all: $(HTML_FILES)

$(HTML_FILES) : %.html : %.qmd
	@echo "Compiling documentation for analyses"
	$(info Source = $<; Destination = $@)
	@echo $(notdir $@)
	$(eval PUBLIC_DIR = $(DOCS_DIR))
	echo "library(quarto); quarto_render(\"$<\")" | R --no-save --no-restore;
	#echo "library(rmarkdown); source(\"resources/render.R\"); rmarkdown::render(\"$<\", output_format = my_html_document(template = \"../resources/tutorial_template.html\", css = \"../resources/tutorial_style.css\", theme = \"spacelab\", toc = TRUE, toc_float = TRUE, highlight = \"pygments\"))" | R --no-save --no-restore;\
	#cp $@ $(DOCS_DIR)$(notdir $@)
