
all: clean dist  install

dist:
	python3 -m build .

clean:
	rm -rf dist *.egg-info

install:
	pip install  --force-reinstall ./dist/*.whl
