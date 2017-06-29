.PHONY: clean-pyc clean-build docs clean updateversion singularity

define BROWSER_PYSCRIPT
import os, webbrowser, sys
try:
        from urllib import pathname2url
except:
        from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT
BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "clean-test - remove test and coverage artifacts"
	@echo "lint - check style with flake8"
	@echo "test - run tests quickly with the default Python"
	@echo "test-all - run tests on every Python version with tox"
	@echo "coverage - check code coverage quickly with the default Python"
	@echo "docs - generate Sphinx HTML documentation, including API docs"
	@echo "release - package and upload a release to pypi"
	@echo "testrelease - package and upload a release to pypitest"
	@echo "dist - package"
	@echo "install - install the package to the active Python's site-packages"
	@echo "updateversion - updates version value in setup.py & d3r/__init__.py"
	@echo "singularity - Creates singularity image"

clean: clean-build clean-pyc clean-test

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -fr {} +

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test:
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/

lint:
	flake8 d3r tests

test:
	python setup.py test

test-all:
	tox

coverage:
	coverage run --source d3r setup.py test
	coverage report -m
	coverage html
	$(BROWSER) htmlcov/index.html

docs:
	rm -f docs/d3r.rst
	rm -f docs/modules.rst
	sphinx-apidoc -o docs/ d3r
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	open docs/_build/html/index.html

testrelease: dist
	twine upload dist/* -r testpypi

release: dist
	twine upload dist/*
	
dist: clean
	python setup.py sdist
	python setup.py bdist_wheel
	ls -l dist

install: clean
	python setup.py install

updateversion:
	@cv=`egrep '^\s+version=' setup.py | sed "s/^.*='//" | sed "s/'.*//"`; \
	read -p "Current ($$cv) enter new version: " vers; \
	echo "Updating setup.py & d3r/__init__.py with new version: $$vers"; \
	sed -i "s/version='.*',/version='$$vers',/" setup.py ; \
	sed -i "s/__version__ = '.*'/__version__ = '$$vers'/" d3r/__init__.py
	@echo -n "  Updated setup.py: " ; \
	grep "version" setup.py ; 
	@echo -n "  Updated d3r/__init__.py: " ; \
	grep "__version__" d3r/__init__.py

singularity: dist
	@echo 'Creating Singularity image'
	vers=`egrep '^\s+version=' setup.py | sed "s/^.*='//" | sed "s/'.*//"`; \
	echo 'version $vers'; \
	imgfile=`echo dist/d3r-$$vers.img` ; \
	whfile=`echo d3r-$$vers-py2.py3-none-any.whl` ; \
	echo 'image file $imgfile' ; \
	sudo singularity create -s 4096 $$imgfile ; \
	sudo singularity bootstrap $$imgfile singularity/d3rcentos.def $$vers; \
	echo 'Singularity Image created $imgfile'
	ls -l dist
