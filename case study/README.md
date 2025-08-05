# Case Study Reproduction Guide

## 1. Data Source

The dataset used in this case study is available at:
[https://amazon-reviews-2023.github.io/data\_processing/0core.html](https://amazon-reviews-2023.github.io/data_processing/0core.html)

We use the **appliances** category. After downloading and extracting the dataset, you will obtain a file named `Appliances.csv`.

### Preprocessing

Place `Appliances.csv` into the `case study` folder. Then run:

```bash
python preprocess_dataset.py
```

This will generate a bipartite graph data file named `amazon_appliances.txt`, which can be read by our program.

---

## 2. Insertion-Only Update Mode

Navigate to the `insertion-only` folder.

### 2.1 Generate Instruction File

Run:

```bash
python generate_instruction_insertion_only.py
```

This will produce `instructions_insonly_appliances_32.txt`, which contains both **insertion** and **query** instructions. By default, each session contains 32 queries (you can modify this value at the beginning of the script `generate_instruction_insertion_only.py`).

### 2.2 Online Algorithm (Online)

Compile and execute:

```bash
g++ ./online_insertion_only.cpp -o online_insertion_only.exe -O3 -std=c++11
./online_insertion_only.exe ../amazon_appliances.txt ./instructions_insonly_appliances_32.txt
```

This runs the Online algorithm on the AM-app dataset for the insertion-only scenario.

### 2.3 Online++ Algorithm

Compile and execute:

```bash
g++ ./onlinepp_insertion_only.cpp -o onlinepp_insertion_only.exe -O3 -std=c++11
./onlinepp_insertion_only.exe ../amazon_appliances.txt ./instructions_insonly_appliances_32.txt
```

This runs the Online++ algorithm on the AM-app dataset for the insertion-only scenario.

### 2.4 BD-Index Algorithm

Compile and execute:

```bash
g++ ./dynamic_insertion_only.cpp -o dynamic_insertion_only.exe -O3 -std=c++11
./dynamic_insertion_only.exe ../amazon_appliances.txt <strategy> ./instructions_insonly_appliances_32.txt
```

Replace `<strategy>` with either `-space` (space-efficient) or `-time` (time-efficient). The turnaround time for queries and updates will be printed to the console.

---

## 3. First-In-First-Out (FIFO) Update Mode

Navigate to the `first-in-first-out` folder.

### 3.1 Generate Subgraph

Since the experiment only inserts **100,000 edges** (much fewer than the original 2.1M edges), first create a subgraph:

```bash
python get_small_graph.py
```

This generates `amazon_appliances_100000.txt`.

### 3.2 Generate Instruction File

Run:

```bash
python generate_instruction_fifo.py
```

This produces `instructions_fifo_256.txt`, which contains both **insertion** and **query** instructions. By default, each session contains 256 queries (you can modify this value at the beginning of the script `generate_instruction_fifo.py`).

### 3.3 Online Algorithm (Online)

Compile and execute:

```bash
g++ ./online_fifo.cpp -o online_fifo.exe -O3 -std=c++11
./online_fifo.exe ./amazon_appliances_100000.txt ./instructions_fifo_256.txt
```

This runs the Online algorithm on the AM-app dataset for the FIFO scenario.

### 3.4 Online++ Algorithm

Compile and execute:

```bash
g++ ./onlinepp_fifo.cpp -o onlinepp_fifo.exe -O3 -std=c++11
./onlinepp_fifo.exe ./amazon_appliances_100000.txt ./instructions_fifo_256.txt
```

This runs the Online++ algorithm on the AM-app dataset for the FIFO scenario.

### 3.5 BD-Index Algorithm

Compile and execute:

```bash
g++ ./dynamic_fifo.cpp -o dynamic_fifo.exe -O3 -std=c++11
./dynamic_fifo.exe ./amazon_appliances_100000.txt <strategy> ./instructions_fifo_256.txt
```

Replace `<strategy>` with either `-space` (space-efficient) or `-time` (time-efficient). The turnaround time for queries and updates will be printed to the console.

---

## 4. Batch-Update and Batch-Query Experiments

For the batch-update and batch-query experiments, simply download the AM dataset from the parent directory and execute the corresponding algorithms provided in this repository.