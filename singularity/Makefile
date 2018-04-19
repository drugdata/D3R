.PHONY: clean singularity

help:
	@echo "clean - remove all build and test artifacts"
	@echo "singularity - Creates singularity 2.3.4 image"

clean:
	rm -fr /home/vagrant/build/


singularity: clean
	@echo 'Creating Singularity v23 image'
	mkdir -p /home/vagrant/build
	imgfile='/home/vagrant/build/d3r.img' ; \
	sudo singularity create -s 12048 $$imgfile ; \
	sudo singularity bootstrap $$imgfile d3rcentos.def; \
	echo 'Singularity image created $$imgfile'
	ls -l /home/vagrant/build
