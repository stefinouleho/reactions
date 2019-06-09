CFLAGS=-g -Wall
CC = gcc -Wall -Wextra
CXX = g++ -Wall -Wextra

run: similarite_cycles
	./similarite_cycles 

init: 
	rm -rf Dossier
	rm -f dernier
	rm -f prochain
	rm -f tmp
	rm -f fichiers/liste_molecules_mces
	rm -f fichiers/liste_molecules_mces.save
	echo "0 0" > prochain
	cp fichiers/molecules.data fichiers/molecules.data.save
	mkdir Dossier
	touch Dossier/similarite.result
	touch Dossier/temps.result

val: similarite_cycles
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./similarite_cycles

similarite_cycles: similarite_cycles.o lecture_molecule_sdf.o graphe_cycles.o mces.o helpers/graph.o helpers/cliquerecursif.o
	gcc ${CFLAGS} -o $@ $^

similarite_cycles_scip: similarite_cycles.o lecture_molecule_sdf.o graphe_cycles.o mces.o helpers/graph.o helpers/proglin_helper_scip.o helpers/sciplib.a helpers/cliquescip.o 
	$(CXX) -I helpers/scip -o $@ $^ -lpopt -lgmp -lm -lz -lreadline -lncurses 	

similarite_cycles.o: similarite_cycles.c similarite_cycles.h structure.h
	gcc ${CFLAGS} -c similarite_cycles.c

graphe_cycles.o: graphe_cycles.c graphe_cycles.h
	gcc ${CFLAGS} -c graphe_cycles.c

mces.o: mces.c mces.h graphe_cycles.h
	gcc ${CFLAGS} -c mces.c		
	
lecture_molecule_sdf.o: lecture_molecule_sdf.c lecture_molecule_sdf.h 
	gcc ${CFLAGS} -c lecture_molecule_sdf.c

helpers/proglin_helper_scip.o: helpers/proglin_helper_scip.c helpers/proglin_helper.h helpers/sciplib.a
	$(CC) -I helpers/scip -o $@ -c $<

clean: 
	rm -f similarite_cycles
	rm -f *o
	rm -f helpers/*.o

