FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors -c
LIBS = -lm
OFILES	=  errorChecking.o symSparseMat.o group.o stack.o network.o clusterDivider.o main.o
EXENAME = cluster

all: $(EXENAME)

clean:
	rm -rf *.o $(EXENAME)

cluster: $(OFILES) 
	gcc $(OFILES) -o $(EXENAME) $(LIBS)

main.o: main.c
	gcc $(FLAGS) main.c

errorChecking.o: errorChecking.c
	gcc $(FLAGS) errorChecking.c
	
clusterDivider.o: clusterDivider.c
	gcc $(FLAGS) clusterDivider.c

group.o: group.c
	gcc $(FLAGS) group.c

network.o: network.c
	gcc $(FLAGS) network.c

stack.o: stack.c
	gcc $(FLAGS) stack.c

symSparseMat.o: symSparseMat.c
	gcc $(FLAGS) symSparseMat.c

