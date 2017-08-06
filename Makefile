.PHONY: all lint test test-cov viz-summarize install dev clean distclean

all: viz-summarize

lint:
	q2lint
	flake8

test: all
	py.test

test-cov: all
	py.test --cov=q2_feature_table

q2_feature_table/_summarize/summarize_assets/dist:
	cd q2_feature_table/_summarize/summarize_assets && \
	npm install && \
	npm run build && \
	cp licenses/* dist

viz-summarize: q2_feature_table/_summarize/summarize_assets/dist

install: all
	python setup.py install

dev: all
	pip install -e .

clean: distclean
	rm -rf q2_feature_table/_summarize/summarize_assets/node_modules

distclean:
	rm -rf q2_feature_table/_summarize/summarize_assets/dist
