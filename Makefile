.PHONY: format lint check all

format:
	isort jbst
	black jbst

lint:
	pylint --exit-zero jbst/seq_tools
	pylint --exit-zero jbst/vector_build



all: format lint
