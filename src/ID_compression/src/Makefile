# Makefile for building C programs

OBJ=sam_models.o sam_file_allocation.o  Arithmetic_stream.o stream_model.o compression.o io_functions.o id_compression.o ../../decompress_quality_id.o #../../reorder_quality_pe.o
CC=g++
RM=rm -f

CXXFLAGS=-O3 -march=native -ffast-math -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -Wall -I../include -I/usr/local/include -I../libssh/build/src  -D__APPLE__ -std=c++11
LDFLAGS=-L/usr/local/lib -lc -lm -lpthread -I../include -std=c++11

%.o : %.c
	$(CC) -gdwarf-2 $(CXXFLAGS) -c $< 
%.o : %.cc
	$(CC) -gdwarf-2 $(CXXFLAGS) -c $< 
%.o : %.cpp
	$(CC) -gdwarf-2 $(CXXFLAGS) -c $< 


all : program

program : $(OBJ)
	$(CC) $(OBJ) -o program $(LDFLAGS)

debug : CXXFLAGS += -DDEBUG -ggdb
debug : program

clean :
	$(RM) *.o program
