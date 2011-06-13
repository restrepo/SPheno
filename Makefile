F90 = ifort
bin/SPheno:
	cd src ; ${MAKE} F90=${F90}
clean:
	rm -f *.o *~ */*.o */*~
cleanall:
	rm -f bin/SPheno lib/*.a *.o *~ */*.o */*~ include/*
.PHONY: bin/SPheno clean cleanall
