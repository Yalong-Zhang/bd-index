import random
from collections import deque

# --------------------------------------------------
# Configuration (modify these parameters as needed)
# --------------------------------------------------
INSERTION_NUMBER = 100000         # Number of edges to insert (use the first N edges from the dataset)
WINDOW_SIZE = 365 * 24 * 60 * 60 * 1000    # Sliding window size in milliseconds (e.g., 1 year)
AB_RANGE = 3                      # Range [0, AB_RANGE] for query parameters a and b
SESSION_NUMBER = 20              # Number of query sessions to simulate
QUERY_NUMBER_PER_SESSION = 256     # Number of queries per session
RANDOM_SEED = 123                 # Random seed for reproducibility

INPUT_FILE = "amazon_appliances_100000.txt"  # Input bipartite dataset file
OUTPUT_FILE = "instructions_fifo_256.txt"        # Output instruction sequence file
SESSION_DURATION_MS = 10 * 1000              # Session duration in milliseconds (10 seconds)
# --------------------------------------------------

def main():
    # 1) Initialize random seed
    random.seed(RANDOM_SEED)

    # 2) Read the bipartite dataset
    with open(INPUT_FILE, 'r') as f:
        # First line: |E| |U| |V|
        header = f.readline().strip().split()
        total_edges, user_count, item_count = map(int, header)

        # Read all edges: u v t
        all_edges = []
        for line in f:
            u_str, v_str, t_str = line.strip().split()
            u, v, t = int(u_str), int(v_str), int(t_str)
            all_edges.append((u, v, t))

    # 3) Take the first INSERTION_NUMBER edges sorted by time
    insertion_edges = sorted(all_edges, key=lambda x: x[2])[:INSERTION_NUMBER]

    # 4) Generate insertion and FIFO deletion instructions
    window = deque()      # holds (u, v, t) within the sliding window
    instructions = []     # will collect all (timestamp, instruction_str)

    for u, v, t in insertion_edges:
        # a) insertion at original timestamp
        instructions.append((t, f"i {u} {v} {t}"))
        # b) add to window
        window.append((u, v, t))
        # c) evict expired edges
        while window and window[0][2] + WINDOW_SIZE <= t:
            u_old, v_old, t_old = window.popleft()
            delete_time = t_old + WINDOW_SIZE
            instructions.append((delete_time, f"d {u_old} {v_old} {delete_time}"))

    # 5) Sort update instructions (insertions + deletions) by timestamp
    instructions.sort(key=lambda x: x[0])

    # 6) Extract only the update instructions for session boundary sampling
    update_insts = [
        (ts, instr) for ts, instr in instructions
        if instr.startswith('i ') or instr.startswith('d ')
    ]

    # 7) Generate query sessions based on adjacent update instructions
    for _ in range(SESSION_NUMBER):
        # pick a random adjacent pair of updates
        idx = random.randint(0, len(update_insts) - 2)
        t_start_window, _ = update_insts[idx]
        t_end_window, _   = update_insts[idx + 1]
        # choose session start uniformly in [t_start_window, t_end_window]
        session_start = random.randint(t_start_window, t_end_window)

        # generate queries within this session
        for _ in range(QUERY_NUMBER_PER_SESSION):
            a = random.randint(0, AB_RANGE)
            b = random.randint(0, AB_RANGE)
            query_time = random.randint(session_start, session_start + SESSION_DURATION_MS)
            instructions.append((query_time, f"q {a} {b} {query_time}"))

    # 8) Final sort of all instructions (updates + queries) by timestamp
    instructions.sort(key=lambda x: x[0])

    # 9) Write out to the instruction file
    with open(OUTPUT_FILE, 'w') as fout:
        fout.write(f"{len(instructions)}\n")
        for _, instr in instructions:
            fout.write(instr + "\n")

    print(f"Generated {len(instructions)} FIFO instructions with interleaved queries in '{OUTPUT_FILE}'")

if __name__ == "__main__":
    main()
