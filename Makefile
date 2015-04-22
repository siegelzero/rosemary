#
# Makefile for rosemary
#

SHELL := /bin/bash
HIDE ?= @
VENV ?= env
PACKAGE := $(shell python setup.py --name)
PIP := $(VENV)/bin/pip

#
# Build rules
#

all:

prepare-venv:
	$(HIDE)virtualenv $(VENV)
	$(HIDE)$(VENV)/bin/pip install --upgrade -r requirements.txt
	$(HIDE)$(VENV)/bin/pip install --upgrade -e .

clean:
	$(HIDE)find . -name \*.pyc -exec rm {} \;


#
# Test rules
#

test: test-unit test-flake8

test-unit:
	$(HIDE)$(VENV)/bin/nosetests $(PACKAGE)

test-flake8:
	$(HIDE)$(VENV)/bin/flake8 --max-line-length=120 $(PACKAGE)

coverage:
	$(HIDE)$(VENV)/bin/nosetests --with-coverage --cover-erase --cover-inclusive --cover-package=$(PACKAGE) $(PACKAGE)
