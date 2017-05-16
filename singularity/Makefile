.PHONY: clean singularity

help:
	@echo "clean - remove all build and test artifacts"
	@echo "singularity - Creates singularity 2.2 image"

clean:
	rm -fr build/


singularity: clean
	@echo 'Creating Singularity v22 image'
	mkdir -p build
	imgfile='build/d3r.img' ; \
	sudo singularity create -s 4096 $$imgfile ; \
	sudo singularity bootstrap $$imgfile d3rcentos.def; \
	echo 'Singularity image created $imgfile'
	ls -l build
