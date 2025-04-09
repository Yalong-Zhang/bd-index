# Optimal (alpha,beta)-Dense Subgraph Search in Static and Dynamic Bipartite Graphs: an Index-Based Approach

# Datasets

Because some datasets used in the paper are too large to be uploaded to GitHub, we have summarized the download links for the dataset in the table below.

Datasets used for performance studies.

| Dataset | Link |
| --- | --- |
| AC | http://www.konect.cc/networks/actor-movie/ |
| IM | http://www.konect.cc/networks/actor2/ |
| HE | http://www.konect.cc/networks/ca-cit-HepPh/ |
| AM | http://www.konect.cc/networks/amazon-ratings/ |
| FL | http://www.konect.cc/networks/flickr-groupmemberships/ |
| EP | http://www.konect.cc/networks/epinions-rating/ |
| PA | http://www.konect.cc/networks/patentcite/ |
| PO | http://www.konect.cc/networks/soc-pokec-relationships/ |
| WI | http://www.konect.cc/networks/edit-eswiki/ |
| LI | http://www.konect.cc/networks/livejournal-groupmemberships/ |

# Preprocess

The dataset needs to be preprocessed into a specific format file to be input by the algorithm, with the format as follows: The first line consists of three integers representing |E|, |U|, |V|. Each subsequent line contains two numbers representing an edge (u, v). Note that |U| should be equal to the maximum value of u among all edges, and |V| should be equal to the maximum value of v among all edges. An example is given in example_graph.txt.

# Usage of Static Algorithm

```
g++ static.cpp -o static -std=c++11 -O3
./static <dataset_address>
```

# Usage of Dynamic Algorithm

```
g++ dynamic.cpp -o dynamic -std=c++11 -O3
./dynamic <dataset_address> <maintenance_strategy>
```

When <maintenance_strategy> is -space, algorithm use the space-efficient maintenance strategy. When <maintenance_strategy> is -time, algorithm use the time-efficient maintenance strategy. 