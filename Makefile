
CXX=g++
CXXFLAGS=-std=gnu++17 -g -fopenmp

default: main

setup: data/804.mxs

data/804.mxs: data/804.mxs.gz
	gunzip -k $<

.PHONY: demo
demo: main
	DEMO=yes ./main data/804.mxs data/Query1.txt demo.json
	cat demo.json

.PHONY:test
test: main  data/804.mxs
	./main data/804.mxs data/Query1.txt 804-Query1.json
	./check.py 804-Query1.json  < data/804-Query1-ref.json
	./main data/804.mxs data/Query5.txt 804-Query5.json
	./check.py 804-Query5.json  < data/804-Query5-ref.json

.PHONY: clean
clean:
	rm -rf main
	rm -rf task.tgz

.PHONY: tgz test-tgz
tgz:task.tgz

task.tgz:
	tar czf $@ $$(git ls-tree -r --name-only master)

test-tgz: task.tgz
	rm -rf task-test
	mkdir task-test
	cd task-test; tar xzf ../$<
	make test
