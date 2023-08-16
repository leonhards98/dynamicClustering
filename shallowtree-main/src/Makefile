INCLUDES = -I boost_1_81_0/boost_1_81_0 -I blaze-3.8.1/
CFLAGS = -c -O3
BuildPath = ../build/

dynC = dynamic_coreset/dyn_coreset.cpp dynamic_coreset/coreset_tree.cpp dynamic_coreset/coreset_node.cpp
cC = coreset/stop_watch.cpp coreset/sensitivity_sampling.cpp coreset/random.cpp coreset/kmeans.cpp coreset/coreset.cpp coreset/clustering_result.cpp coreset/cluster_assignment_list.cpp
coC = coreset_old/stop_watch.cpp coreset_old/sensitivity_sampling.cpp coreset_old/random.cpp coreset_old/kmeans.cpp coreset_old/coreset.cpp coreset_old/clustering_result.cpp coreset_old/cluster_assignment_list.cpp
kmeansC = kmeanspp/fastKMeans.cpp kmeanspp/sampleTree.cpp kmeanspp/quadTree.cpp kmeanspp/quadTreeNode.cpp
ILPC = ILPBaseLine/baseline.cpp

dynObj = $(BuildPath)dynamic_coreset/dyn_coreset.o $(BuildPath)dynamic_coreset/coreset_tree.o $(BuildPath)dynamic_coreset/coreset_node.o
cObj = $(BuildPath)coreset/stop_watch.o $(BuildPath)coreset/sensitivity_sampling.o $(BuildPath)coreset/random.o $(BuildPath)coreset/kmeans.o $(BuildPath)coreset/coreset.o $(BuildPath)coreset/clustering_result.o $(BuildPath)coreset/cluster_assignment_list.o
coObj = $(BuildPath)coreset_old/stop_watch.o $(BuildPath)coreset_old/sensitivity_sampling.o $(BuildPath)coreset_old/random.o $(BuildPath)coreset_old/kmeans.o $(BuildPath)coreset_old/coreset.o $(BuildPath)coreset_old/clustering_result.o $(BuildPath)coreset_old/cluster_assignment_list.o
kmeansObj = $(BuildPath)kmeanspp/fastKMeans.o $(BuildPath)kmeanspp/sampleTree.o $(BuildPath)kmeanspp/quadTree.o $(BuildPath)kmeanspp/quadTreeNode.o
ILPObj = $(BuildPath)ILPBaseLine/baseline.o

#$(dynObj): $(BUILDDIR)/%.o : $(dynC)
	#$(CC) $(CFLAGS) $(INCLUDES) $< -o $@

main: $(dynObj) $(cObj) $(coObj) $(kmeansObj)
	g++ main.cpp $(dynObj) $(cObj) $(coObj) $(kmeansObj) $(INCLUDES) -o $(BuildPath)a.out



$(dynObj): $(dynC)
	g++ -c $(patsubst $(BuildPath)%.o,%.cpp,$@) $(INCLUDES) $(FLAGS) -o $@

$(cObj): $(cC)
	g++ -c $(patsubst $(BuildPath)%.o,%.cpp,$@) $(INCLUDES) $(FLAGS) -o $@

$(coObj): $(coC)
	g++ -c $(patsubst $(BuildPath)%.o,%.cpp,$@) $(INCLUDES) $(FLAGS) -o $@

$(kmeansObj): $(kmeansC)
	g++ -c $(patsubst $(BuildPath)%.o,%.cpp,$@) $(INCLUDES) $(FLAGS) -o $@

clear:
	rm $(BuildPath)dynamic_coreset/*.o
	rm $(BuildPath)coreset/*.o
	rm $(BuildPath)coreset_old/*.o
	rm $(BuildPath)kmeanspp/*.o
	rm a.out

