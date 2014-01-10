#
# CDDL HEADER START
#
# The contents of this file are subject to the terms of the Common Development
# and Distribution License Version 1.0 (the "License").
#
# You can obtain a copy of the license at
# http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
# specific language governing permissions and limitations under the License.
#
# When distributing Covered Code, include this CDDL HEADER in each file and
# include the License file in a prominent location with the name LICENSE.CDDL.
# If applicable, add the following below this CDDL HEADER, with the fields
# enclosed by brackets "[]" replaced with your own identifying information:
#
# Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
#
# CDDL HEADER END
#

#
# Copyright (c) 2013, Regents of the University of Minnesota.
# All rights reserved.
#
# Contributors:
#    Ryan S. Elliott
#    Valeriu Smirichinski
#    Ellad B. Tadmor
#

#
# Release: This file is part of the openkim-api.git repository.
#


include Makefile.KIM_Config

KIM_LISTS = Makefile.ModelsList Makefile.ModelDriversList Makefile.TestsList

KIM_CONFIG_FILES = $(KIM_DIR)/KIM_API/Makefile.KIM_Config $(KIM_MODEL_DRIVERS_DIR)/Makefile.KIM_Config $(KIM_MODELS_DIR)/Makefile.KIM_Config $(KIM_TESTS_DIR)/Makefile.KIM_Config

include $(KIM_LISTS)

.PHONY: all models_check config kim_lists \
            $(patsubst %,%-all,  $(MODELS_LIST) $(MODEL_DRIVERS_LIST) $(TESTS_LIST)) \
        clean kim-api-clean config-clean kim_lists_clean \
            $(patsubst %,%-clean,$(MODELS_LIST) $(MODEL_DRIVERS_LIST) $(TESTS_LIST)) \
        install kim-api-install config-install \
            $(patsubst %,%-install,$(MODELS_LIST) $(MODEL_DRIVERS_LIST)) \
        uninstall kim-api-uninstall config-uninstall \
            $(patsubst %,%-uninstall,$(MODELS_LIST) $(MODEL_DRIVERS_LIST)) \
        openkim-api kim-api-objects kim-api-libs \
        examples examples-all \
        examples-force examples-force-all \
        examples-clean examples-clean-all

# compile everything in the standard directories
ifeq ($(KIM_LINK),dynamic-load)
   all: models_check config kim-api-objects kim-api-libs $(patsubst %,%-all,$(MODEL_DRIVERS_LIST) \
        $(MODELS_LIST)) $(patsubst %,%-all,$(TESTS_LIST))
else # everything else: dynamic-link static-link
   all: models_check config kim-api-objects $(patsubst %,%-all,$(MODEL_DRIVERS_LIST) $(MODELS_LIST)) \
        kim-api-libs $(patsubst %,%-all,$(TESTS_LIST))
endif

# other targets
clean: config $(patsubst %,%-clean,$(MODELS_LIST) $(MODEL_DRIVERS_LIST) $(TESTS_LIST)) kim-api-clean kim_lists-clean config-clean
install: config $(patsubst %,%-install,$(MODELS_LIST) $(MODEL_DRIVERS_LIST)) kim-api-install config-install
uninstall: config $(patsubst %,%-uninstall,$(MODELS_LIST) $(MODEL_DRIVERS_LIST)) kim-api-uninstall config-uninstall uninstall-cleanup
openkim-api: config kim-api-objects kim-api-libs # compile the openkim-api
examples: config examples-all                    # copy examples to appropriate directories then make
examples-force: config examples-force-all
examples-clean: examples-clean-all

########### for internal use ###########
%-making-echo:
	@printf "\n%79s\n" " " | sed -e 's/ /*/g'
	@printf "%-77s%2s\n" "** Making... `printf "$(patsubst %-all,%,$*)" | sed -e 's/@/ /g'`" "**"
	@printf "%79s\n" " " | sed -e 's/ /*/g'

config: $(KIM_CONFIG_FILES)

kim_lists: $(KIM_LISTS)

kim_lists-clean:
	@printf "Cleaning... KIM_LISTS files.\n"
	@rm -f $(KIM_LISTS)

Makefile: $(KIM_LISTS)
	@touch Makefile

Makefile.ModelsList: $(KIM_MODELS_DIR)
	@printf "Creating... $@.\n"
	@printf "MODELS_LIST = $(notdir $(filter-out $(shell if test -e "$(KIM_MODELS_DIR)/.kimignore"; then cat "$(KIM_MODELS_DIR)/.kimignore";fi;),$(filter-out .%,$(shell find "$(KIM_MODELS_DIR)/" -maxdepth 1 -mindepth 1 -type d -exec basename {} \;))))\n" > Makefile.ModelsList

Makefile.ModelDriversList: $(KIM_MODEL_DRIVERS_DIR)
	@printf "Creating... $@.\n"
	@printf "MODEL_DRIVERS_LIST = $(notdir $(filter-out $(shell if test -e "$(KIM_MODEL_DRIVERS_DIR)/.kimignore"; then cat "$(KIM_MODEL_DRIVERS_DIR)/.kimignore";fi;),$(filter-out .%,$(shell find "$(KIM_MODEL_DRIVERS_DIR)/" -maxdepth 1 -mindepth 1 -type d -exec basename {} \;))))\n" > Makefile.ModelDriversList

Makefile.TestsList: $(KIM_TESTS_DIR)
	@printf "Creating... $@.\n"
	@printf "TESTS_LIST = $(notdir $(filter-out $(shell if test -e "$(KIM_TESTS_DIR)/.kimignore"; then cat "$(KIM_TESTS_DIR)/.kimignore";fi;),$(filter-out .%,$(shell find "$(KIM_TESTS_DIR)/" -maxdepth 1 -mindepth 1 -type d -exec basename {} \;))))\n" > Makefile.TestsList

$(KIM_CONFIG_FILES): Makefile $(KIM_DIR)/Makefile.KIM_Config
	@printf "Creating... KIM_Config file..... $(patsubst $(KIM_DIR)/%,%,$@).\n";
	@printf "# This file is automatically generated by the openkim-api make system.\n" > $@; \
         printf "# Do not edit!\n"                                                        >> $@; \
         printf "\n"                                                                      >> $@; \
         printf "include $(KIM_DIR)/Makefile.KIM_Config\n"                                 >> $@;

config-clean:
	@printf "Cleaning... KIM_Config files.\n"
	@rm -f $(KIM_CONFIG_FILES)

config-install:
	@printf "Installing... KIM_Config files"
ifeq (dynamic-load,$(KIM_LINK))
	@printf ".\n"
	@cp -rf MAKE_SYSTEM $(DESTDIR)$(libdir)/$(package_name)/MAKE_SYSTEM
	@chmod 755 $(DESTDIR)$(libdir)/$(package_name)/MAKE_SYSTEM
	@chmod 755 $(DESTDIR)$(libdir)/$(package_name)/MAKE_SYSTEM/*_DEFAULTS
	@chmod 644 $(DESTDIR)$(libdir)/$(package_name)/MAKE_SYSTEM/*_DEFAULTS/*
	@chmod 644 $(DESTDIR)$(libdir)/$(package_name)/MAKE_SYSTEM/Makefile.*
	@chmod 644 $(DESTDIR)$(libdir)/$(package_name)/MAKE_SYSTEM/model_based_on_driver.cpp
	@sed -e 's|^ *prefix *=.*|prefix = $(prefix)|' \
             -e 's|^ *libdir *=.*|libdir = $(libdir)|' \
             -e '/ *#INSTALLTAG/,/ *#INSTALLTAG/d'     \
             MAKE_SYSTEM/Makefile.Generic > $(DESTDIR)$(libdir)/$(package_name)/MAKE_SYSTEM/Makefile.Generic
	@chmod 644 $(DESTDIR)$(libdir)/$(package_name)/MAKE_SYSTEM/Makefile.Generic
	@sed -e '/KIM_MODEL_DRIVERS_DIR/d' \
             -e '/KIM_MODELS_DIR/d'        \
             -e '/KIM_TESTS_DIR/d'         \
             -e 's|^ *KIM_DIR *=.*|KIM_DIR = $(libdir)/$(package_name)|' \
             Makefile.KIM_Config > $(DESTDIR)$(libdir)/$(package_name)/Makefile.KIM_Config
	@chmod 644 $(DESTDIR)$(libdir)/$(package_name)/Makefile.KIM_Config
else ifeq (dynamic-link,$(KIM_LINK))
	@if test \( -d $(DESTDIR)$(libdir)/$(package_name)/MAKE_SYSTEM \) -o \( -f $(DESTDIR)$(libdir)/$(package_name)/Makefile.KIM_Config \); then \
            printf ": removing KIM_Config files for dynamic-load"; \
            rm -rf "$(DESTDIR)$(libdir)/$(package_name)/MAKE_SYSTEM"; \
            rm -f "$(DESTDIR)$(libdir)/$(package_name)/Makefile.KIM_Config"; \
         else \
            printf ": nothing to be done for dynamic-link"; \
         fi;
	@printf ".\n"
else
	@printf ": nothing to be done for static-link.\n"
endif

config-uninstall:
	@printf "Uninstalling... KIM_Config files"
ifeq (dynamic-load,$(KIM_LINK))
	@printf ".\n"
	@if test \( -d "$(DESTDIR)$(libdir)/$(package_name)/MAKE_SYSTEM" \); then rm -rf "$(DESTDIR)$(libdir)/$(package_name)/MAKE_SYSTEM"; fi;
	@if test \( -f "$(DESTDIR)$(libdir)/$(package_name)/Makefile.KIM_Config" \); then rm -f "$(DESTDIR)$(libdir)/$(package_name)/Makefile.KIM_Config"; fi;
else ifeq (dynamic-link,$(KIM_LINK))
	@printf ": nothing to be done for dynamic-link.\n"
else
	@printf ": nothing to be done for static-link.\n"
endif


uninstall-cleanup:
	@printf "Uninstalling... $(package_name) directories.\n"
	@if test \( -d "$(DESTDIR)$(libdir)/$(package_name)/MODELS" \); then rmdir "$(DESTDIR)$(libdir)/$(package_name)/MODELS" >& /dev/null; fi
	@if test \( -d "$(DESTDIR)$(libdir)/$(package_name)/MODEL_DRIVERS" \); then rmdir "$(DESTDIR)$(libdir)/$(package_name)/MODEL_DRIVERS" >& /dev/null; fi
	@if test \( -d "$(DESTDIR)$(libdir)/$(package_name)" \); then rmdir "$(DESTDIR)$(libdir)/$(package_name)" >& /dev/null; fi

kim-api-objects: Makefile kim-api-objects-making-echo
	@$(MAKE) $(MAKE_FLAGS) -C $(KIM_DIR)/KIM_API/ objects

kim-api-libs: Makefile kim-api-libs-making-echo
	@$(MAKE) $(MAKE_FLAGS) -C $(KIM_DIR)/KIM_API/ libs

kim-api-clean:
	@$(MAKE) $(MAKE_FLAGS) -C $(KIM_DIR)/KIM_API/ clean
	@rm -f kim.log

kim-api-install:
	@$(MAKE) $(MAKE_FLAGS) -C $(KIM_DIR)/KIM_API/ install

kim-api-uninstall:
	@$(MAKE) $(MAKE_FLAGS) -C $(KIM_DIR)/KIM_API/ uninstall

examples-all:
	@$(foreach exmpl,$(notdir $(shell find $(KIM_DIR)/EXAMPLES/MODEL_DRIVERS -maxdepth 1 -mindepth 1 \( -type d -o -type f \) -exec basename {} \;)),\
          if test -e $(KIM_MODEL_DRIVERS_DIR)/$(exmpl); then \
          printf "*@existing.....@%-50s@no@copy@performed!\n" $(exmpl)@ | sed -e 's/ /./g' -e 's/@/ /g'; else \
          printf "*@installing...@%-50s@copied@to@$(KIM_MODEL_DRIVERS_DIR)\n" $(exmpl)@ | sed -e 's/ /./g' -e 's/@/ /g'; \
          cp -r $(KIM_DIR)/EXAMPLES/MODEL_DRIVERS/$(exmpl) "$(KIM_MODEL_DRIVERS_DIR)/"; fi;)
	@$(foreach exmpl,$(notdir $(shell find $(KIM_DIR)/EXAMPLES/MODELS -maxdepth 1 -mindepth 1 \( -type d -o -type f \) -exec basename {} \;)),\
          if test -e $(KIM_MODELS_DIR)/$(exmpl); then \
          printf "*@existing.....@%-50s@no@copy@performed!\n" $(exmpl)@ | sed -e 's/ /./g' -e 's/@/ /g'; else \
          printf "*@installing...@%-50s@copied@to@$(KIM_MODELS_DIR)\n" $(exmpl)@ | sed -e 's/ /./g' -e 's/@/ /g'; \
          cp -r $(KIM_DIR)/EXAMPLES/MODELS/$(exmpl) "$(KIM_MODELS_DIR)/"; fi;)
	@$(foreach exmpl,$(notdir $(shell find $(KIM_DIR)/EXAMPLES/TESTS -maxdepth 1 -mindepth 1 \( -type d -o -type f \) -exec basename {} \;)),\
          if test -e $(KIM_TESTS_DIR)/$(exmpl); then \
          printf "*@existing.....@%-50s@no@copy@performed!\n" $(exmpl)@ | sed -e 's/ /./g' -e 's/@/ /g'; else \
          printf "*@installing...@%-50s@copied@to@$(KIM_TESTS_DIR)\n" $(exmpl)@ | sed -e 's/ /./g' -e 's/@/ /g'; \
          cp -r $(KIM_DIR)/EXAMPLES/TESTS/$(exmpl) "$(KIM_TESTS_DIR)/"; fi;)

examples-clean: kim_lists-clean
	@printf "Removing all installed example files...\n"
	@$(foreach dr,$(notdir $(wildcard $(KIM_DIR)/EXAMPLES/MODEL_DRIVERS/*)), rm -rf "$(KIM_MODEL_DRIVERS_DIR)/$(dr)";)
	@$(foreach dr,$(notdir $(wildcard $(KIM_DIR)/EXAMPLES/MODELS/*)), rm -rf "$(KIM_MODELS_DIR)/$(dr)";)
	@$(foreach dr,$(notdir $(wildcard $(KIM_DIR)/EXAMPLES/TESTS/*)), rm -rf "$(KIM_TESTS_DIR)/$(dr)";)

examples-force-all:
	@$(foreach exmpl,$(notdir $(shell find $(KIM_DIR)/EXAMPLES/MODEL_DRIVERS -maxdepth 1 -mindepth 1 \( -type d -o -type f \) -exec basename {} \;)),\
          if test -e $(KIM_MODEL_DRIVERS_DIR)/$(exmpl); then \
          printf "*@overwriting..@%-50s@copied@to@$(KIM_MODEL_DRIVERS_DIR)\n" $(exmpl)@ | sed -e 's/ /./g' -e 's/@/ /g'; \
          rm -rf "$(KIM_MODEL_DRIVERS_DIR)/$(exmpl)"; \
          cp -r $(KIM_DIR)/EXAMPLES/MODEL_DRIVERS/$(exmpl) "$(KIM_MODEL_DRIVERS_DIR)/"; else \
          printf "*@installing...@%-50s@copied@to@$(KIM_MODEL_DRIVERS_DIR)\n" $(exmpl)@ | sed -e 's/ /./g' -e 's/@/ /g'; \
          cp -r $(KIM_DIR)/EXAMPLES/MODEL_DRIVERS/$(exmpl) "$(KIM_MODEL_DRIVERS_DIR)/"; fi;)
	@$(foreach exmpl,$(notdir $(shell find $(KIM_DIR)/EXAMPLES/MODELS -maxdepth 1 -mindepth 1 \( -type d -o -type f \) -exec basename {} \;)),\
          if test -e $(KIM_MODELS_DIR)/$(exmpl); then \
          printf "*@overwriting..@%-50s@copied@to@$(KIM_MODELS_DIR)\n" $(exmpl)@ | sed -e 's/ /./g' -e 's/@/ /g'; \
          rm -rf "$(KIM_MODELS_DIR)/$(exmpl)"; \
          cp -r $(KIM_DIR)/EXAMPLES/MODELS/$(exmpl) "$(KIM_MODELS_DIR)/"; else \
          printf "*@installing..@%-50s@copied@to@$(KIM_MODELS_DIR)\n" $(exmpl)@ | sed -e 's/ /./g' -e 's/@/ /g'; \
          cp -r $(KIM_DIR)/EXAMPLES/MODELS/$(exmpl) "$(KIM_MODELS_DIR)/"; fi;)
	@$(foreach exmpl,$(notdir $(shell find $(KIM_DIR)/EXAMPLES/TESTS -maxdepth 1 -mindepth 1 \( -type d -o -type f \) -exec basename {} \;)),\
          if test -e $(KIM_TESTS_DIR)/$(exmpl); then \
          printf "*@overwriting..@%-50s@copied@to@$(KIM_TESTS_DIR)\n" $(exmpl)@ | sed -e 's/ /./g' -e 's/@/ /g'; \
          rm -rf "$(KIM_TESTS_DIR)/$(exmpl)"; \
          cp -r $(KIM_DIR)/EXAMPLES/TESTS/$(exmpl) "$(KIM_TESTS_DIR)/"; else \
          printf "*@installing..@%-50s@copied@to@$(KIM_TESTS_DIR)\n" $(exmpl)@ | sed -e 's/ /./g' -e 's/@/ /g'; \
          cp -r $(KIM_DIR)/EXAMPLES/TESTS/$(exmpl) "$(KIM_TESTS_DIR)/"; fi;)

$(patsubst %,%-all,$(MODELS_LIST)): %: Makefile Model..........@%-making-echo | kim-api-objects
	@$(MAKE) $(MAKE_FLAGS) -C $(KIM_MODELS_DIR)/$(patsubst %-all,%,$@) all

$(patsubst %,%-clean,$(MODELS_LIST)):
	@$(MAKE) $(MAKE_FLAGS) -C $(KIM_MODELS_DIR)/$(patsubst %-clean,%,$@) clean

$(patsubst %,%-install,$(MODELS_LIST)):
	@$(MAKE) $(MAKE_FLAGS) -C $(KIM_MODELS_DIR)/$(patsubst %-install,%,$@) install

$(patsubst %,%-uninstall,$(MODELS_LIST)):
	@$(MAKE) $(MAKE_FLAGS) -C $(KIM_MODELS_DIR)/$(patsubst %-uninstall,%,$@) uninstall

$(patsubst %,%-all,$(MODEL_DRIVERS_LIST)): %: Makefile Model@Driver...@%-making-echo | kim-api-objects
	@$(MAKE) $(MAKE_FLAGS) -C $(KIM_MODEL_DRIVERS_DIR)/$(patsubst %-all,%,$@) all

$(patsubst %,%-clean,$(MODEL_DRIVERS_LIST)):
	@$(MAKE) $(MAKE_FLAGS) -C $(KIM_MODEL_DRIVERS_DIR)/$(patsubst %-clean,%,$@) clean

$(patsubst %,%-install,$(MODEL_DRIVERS_LIST)):
	@$(MAKE) $(MAKE_FLAGS) -C $(KIM_MODEL_DRIVERS_DIR)/$(patsubst %-install,%,$@) install

$(patsubst %,%-uninstall,$(MODEL_DRIVERS_LIST)):
	@$(MAKE) $(MAKE_FLAGS) -C $(KIM_MODEL_DRIVERS_DIR)/$(patsubst %-uninstall,%,$@) uninstall

$(patsubst %,%-all,$(TESTS_LIST)): %: Makefile Test...........@%-making-echo | kim-api-objects kim-api-libs $(patsubst %,%-all,$(MODELS_LIST))
	@$(MAKE) $(MAKE_FLAGS) -C $(KIM_TESTS_DIR)/$(patsubst %-all,%,$@) all

$(patsubst %,%-clean,$(TESTS_LIST)):
	@$(MAKE) $(MAKE_FLAGS) -C $(KIM_TESTS_DIR)/$(patsubst %-clean,%,$@) clean

models_check:
	@if test \( X"$(MODELS_LIST)" = X"" \) -a \( X"$(KIM_LINK)" = X"static-link" \); then     \
        printf "*******************************************************************************\n"; \
        printf "*******     Can't compile the API for static linking with no Models     *******\n"; \
        printf "*******              Maybe you want to do 'make examples'               *******\n"; \
        printf "*******************************************************************************\n"; \
        false; else true; fi
