# generate_new_bipartite.py

# --------------------------------------------------
# Configuration
# --------------------------------------------------
INSERTION_NUMBER = 100000        # Number of edges to keep from the original file
INPUT_FILE = "../amazon_appliances.txt"
OUTPUT_FILE = "amazon_appliances_100000.txt"
# --------------------------------------------------

def main():
    # 1) Read the original bipartite file
    with open(INPUT_FILE, 'r') as fin:
        # Read header: total_edges, user_count, item_count
        header = fin.readline().strip().split()
        total_edges, orig_user_count, orig_item_count = map(int, header)
        fin.readline().strip()

        # 2) Read only the first INSERTION_NUMBER edges
        edges = []
        for _ in range(min(INSERTION_NUMBER, total_edges)):
            line = fin.readline().strip()
            if not line:
                break
            u_str, v_str, t_str = line.split()
            u, v, t = int(u_str), int(v_str), int(t_str)
            edges.append((u, v, t))

    # 3) Build new mappings for users and items (separate reindexing)
    user_set = {u for (u, _, _) in edges}
    item_set = {v for (_, v, _) in edges}

    user_list = sorted(user_set)
    item_list = sorted(item_set)

    user2new = {old_u: new_u + 1 for new_u, old_u in enumerate(user_list)}
    item2new = {old_v: new_v + 1 for new_v, old_v in enumerate(item_list)}

    # 4) Write out the new bipartite graph file
    with open(OUTPUT_FILE, 'w') as fout:
        new_edge_count = len(edges)
        new_user_count = len(user_list)
        new_item_count = len(item_list)
        # Write header: |E'| |U'| |V'|
        fout.write(f"{new_edge_count} {new_user_count} {new_item_count}\n")
        # Write each edge with reindexed nodes
        for u_old, v_old, t in edges:
            u_new = user2new[u_old]
            v_new = item2new[v_old]
            fout.write(f"{u_new} {v_new} {t}\n")

    print(f"Wrote {new_edge_count} edges to '{OUTPUT_FILE}'")
    print(f"Users reindexed: {new_user_count}, Items reindexed: {new_item_count}")

if __name__ == "__main__":
    main()
