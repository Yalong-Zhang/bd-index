import random

# --------------------------------------------------
# Configuration (modify these parameters as needed)
# --------------------------------------------------
INSERTION_EDGE_NUMBER = 2000       # Number of edges to insert (use the last N edges from the dataset)
AB_RANGE = 3                       # Range [0, AB_RANGE] for query parameters a and b
SESSION_NUMBER = 20                # Number of query sessions to simulate
QUERY_NUMBER_PER_SESSION = 32      # Number of queries per session
RANDOM_SEED = 123                   # Random seed for reproducibility

INPUT_FILE = "../amazon_appliances.txt"  # Input bipartite dataset file
OUTPUT_FILE = "instructions_insonly_appliances_32.txt"      # Output instruction sequence file
SESSION_DURATION_MS = 10 * 1000       # Session duration in milliseconds (10 seconds)
# --------------------------------------------------

def main():
    # Initialize random seed
    random.seed(RANDOM_SEED)

    # 1) Read the bipartite dataset
    with open(INPUT_FILE, 'r') as f:
        # First line: |E| |U| |V|
        header = f.readline().strip().split()
        total_edges, user_count, item_count = map(int, header)

        # Subsequent lines: u v t
        edges = []
        for line in f:
            u_str, v_str, t_str = line.strip().split()
            u, v, t = int(u_str), int(v_str), int(t_str)
            edges.append((u, v, t))

    # 2) Select the last INSERTION_EDGE_NUMBER edges for insertion commands
    insertion_edges = edges[-INSERTION_EDGE_NUMBER:]
    instructions = []
    for u, v, t in insertion_edges:
        # Format: "i u v t" means insert edge (u,v) at time t
        instructions.append((t, f"i {u} {v} {t}"))

    # 3) Determine the time range of inserted edges
    times = [t for (_, _, t) in insertion_edges]
    min_time, max_time = min(times), max(times)

    # 4) Generate query sessions
    #    Each session has a random start time in [min_time, max_time],
    #    lasts SESSION_DURATION_MS, and contains QUERY_NUMBER_PER_SESSION queries.
    session_starts = [
        random.randint(min_time, max_time)
        for _ in range(SESSION_NUMBER)
    ]
    for start in session_starts:
        for _ in range(QUERY_NUMBER_PER_SESSION):
            # Randomly pick query parameters a and b in [0, AB_RANGE]
            a = random.randint(0, AB_RANGE)
            b = random.randint(0, AB_RANGE)
            # Randomly pick a timestamp within the session window
            query_time = random.randint(start, start + SESSION_DURATION_MS)
            # Format: "q a b t" means query (a,b) at time t
            instructions.append((query_time, f"q {a} {b} {query_time}"))

    # 5) Sort all instructions by timestamp
    instructions.sort(key=lambda x: x[0])

    # 6) Write instructions to the output file
    with open(OUTPUT_FILE, 'w') as fout:
        # First line: total number of instructions
        fout.write(f"{len(instructions)}\n")
        # Subsequent lines: one instruction per line
        for _, instr in instructions:
            fout.write(instr + "\n")

    print(f"Generated {len(instructions)} instructions and saved to '{OUTPUT_FILE}'")

if __name__ == "__main__":
    main()
