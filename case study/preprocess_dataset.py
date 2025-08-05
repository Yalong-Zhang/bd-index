import pandas as pd

def convert_amazon_to_bipartite(input_file, output_file):
    df = pd.read_csv(input_file)

    user2id = {uid: idx + 1 for idx, uid in enumerate(df['user_id'].unique())}
    item2id = {iid: idx + 1 for idx, iid in enumerate(df['parent_asin'].unique())}

    df['user_int'] = df['user_id'].map(user2id).astype(int)
    df['item_int'] = df['parent_asin'].map(item2id).astype(int)

    df = df.sort_values(by='timestamp')

    edge_count = len(df)
    user_count = len(user2id)
    item_count = len(item2id)

    with open(output_file, 'w') as f:
        f.write(f"{edge_count} {user_count} {item_count}\n")
        for _, row in df.iterrows():
            f.write(f"{row['user_int']} {row['item_int']} {int(row['timestamp'])}\n")

if __name__ == "__main__":
    input_file = "Appliances.csv"
    output_file = "amazon_appliances.txt"
    convert_amazon_to_bipartite(input_file, output_file)
