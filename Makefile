SUBDIRS = ./src/bo_python/src \
          ./src/DRV_prediction_model/feature_extractor_c \
          ./src/DRV_prediction_model/innovus_feature_extractor_c

# Default target: run cmake and make in all subdirectories
all: $(SUBDIRS)

# Define a target for each subdirectory
$(SUBDIRS):
	@echo "Building in $@"
	cd $@ && cmake . && $(MAKE)

# Clean target: run 'make clean' in all subdirectories
clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done

.PHONY: all clean $(SUBDIRS)

