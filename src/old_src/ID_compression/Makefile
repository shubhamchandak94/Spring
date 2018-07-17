all:
	echo $(MAKE)
	$(MAKE) -C src
	mkdir -p bin
	mv src/program bin/program

debug:
	$(MAKE) -C -g src debug
	mkdir -p bin
	mv src/program bin/program

clean:
	$(MAKE) -C src clean
	rm -f bin/program
