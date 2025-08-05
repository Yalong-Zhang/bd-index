#include <bits/stdc++.h>
#include <unordered_set>
#include <unordered_map>
#include <chrono>
using namespace std;
const int INF = 2000000000;
typedef long long ll;

inline int read_number(FILE* in) {
	int x = 0; char ch = 0; while (ch < '0' || ch > '9') ch = fgetc(in); while (ch >= '0' && ch <= '9') { x = x * 10 + (ch - '0'); ch = fgetc(in); } return x;
}
inline void check(bool flag, const char* message) {
	if (!flag) {
		printf("!!!!! CHECK ERROR !!!!!\n");
		printf("Error message: %s\n", message);
		assert(0);
	}
}
enum Algorithm_used {
	SPACE, TIME
};
Algorithm_used algorithm_used;

struct Timer {
	chrono::high_resolution_clock::time_point start_time, end_time;
	void start() { start_time = chrono::high_resolution_clock::now(); }
	void end() { end_time = chrono::high_resolution_clock::now(); }
	double time() {
		return chrono::duration<double>(end_time - start_time).count();
	}
};
Timer timer;
template <class T>
struct Set {
	T* nodes; bool* in; int size = -1;
	Set() {}
	Set(int sz) { size = 0; nodes = (T*)malloc(sz * sizeof(T)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); }
	void alloc(int sz) { size = 0; nodes = (T*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); }
	void insert(T x) { nodes[size++] = x; in[x] = true; }
	void clear() { for (int i = 0; i < size; i++) in[nodes[i]] = false; size = 0; }
	~Set() { free(nodes), free(in); }
};
template <class T>
struct Map {
	int* nodes; bool* in; int size = -1; T* value;
	Map() {}
	Map(int sz) { size = 0; nodes = (int*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); value = (T*)malloc(sz * sizeof(T)); memset(value, 0, sz * sizeof(T)); }
	void alloc(int sz) { size = 0; nodes = (int*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); value = (T*)malloc(sz * sizeof(T)); memset(value, 0, sz * sizeof(T)); }
	void freememory() { free(nodes), free(in), free(value); }
	void clear() { for (int i = 0; i < size; i++) in[nodes[i]] = false, value[nodes[i]] = 0; size = 0; }
	T& operator[](int x) { if (!in[x]) nodes[size++] = x, in[x] = true; return value[x]; }
	~Map() { free(nodes), free(in), free(value); }
};
template <class T>
struct Queue {
	T* nodes; int head, tail;
	Queue() {}
	Queue(int size) { nodes = (T*)malloc(size * sizeof(T)); head = tail = 0; }
	void alloc(int sz) { head = tail = 0; nodes = (int*)malloc(sz * sizeof(int)); }
	~Queue() { free(nodes); }
	bool empty() { return head == tail; }
	int pop() { return nodes[head++]; }
	void push(T x) { nodes[tail++] = x; }
	void clear() { head = tail = 0; }
};

struct Edge { int u, v, to; bool deleted; };
struct Graph {
	int M, U, V, N;
	Edge* e;
	int* undeg, * indeg;
	int** adj;
	int pseudoarboricity;
	void read_graph_from_dataset(char* dataset_address);

	inline bool in_U(int x) { return x < U; }
	inline bool in_V(int x) { return x >= U; }
	inline bool in_S(int x) { return indeg[x] < (in_U(x) ? alpha : beta); }
	inline bool in_T(int x) { return indeg[x] > (in_U(x) ? alpha : beta); }

	void initialize_orientation(int ab, bool the_first_iteration);
	int alpha, beta;

	Set<int> C1, C2;
	Set<int> C1_minus_C2;
	void get_core(int core_alpha, int core_beta, Set<int>& C);

	int* sorted; Set<int> D;
	bool* computed; int* position, * number_of_edges, * number_of_nodes;
	void change_sorted(int position_l, int position_u, Set<int>& D);

	int alpha_max, beta_max, d_max;
	void construct_index();
	void ReTest(int _alpha, int _beta, bool the_first_iteration);
	void Divide(int D_l_rank, int D_u_rank, bool the_first_iteration);

	int* r;

	Set<int> S; Queue<int> Q;
	Map<int> parent; Set<int> vis; Map<int> dist, cur;

	bool DinicBFS(int D_l_rank, int D_u_rank);
	bool DinicDFS(int x, int D_l_rank, int D_u_rank);

	void analyze_r();
	void check_correctness();
	const int OUTPUT_NODES_LIMIT = 10;
	void output_core(Set<int>& C);
	void output_dense_subgraph();
	void display();

	void Dynamic_ReTest(int _alpha, int _beta);
	bool Dynamic_DinicBFS();
	bool Dynamic_DinicDFS(int x);

	void delete_edge_space(int edge_id);
	void insert_edge_space(int edge_id);

	void delete_edge_time(int edge_id);
	void insert_edge_time(int edge_id);
};
Graph G;

struct BD_Index {
	vector<vector<list<int>::iterator> > pointer_U;
	vector<list<int> > index_U;
	vector<vector<list<int>::iterator> > pointer_V;
	vector<list<int> > index_V;
	void compute_memory();
	void inc_U(int alpha, int beta, unordered_set<int>& update_nodes);
	void inc_V(int beta, int alpha, unordered_set<int>& update_nodes);
	void dec_U(int alpha, int beta, unordered_set<int>& update_nodes);
	void dec_V(int beta, int alpha, unordered_set<int>& update_nodes);
	void set_U(int alpha, int x);
	void set_V(int beta, int x);
	void set_U_time(int alpha, int x);
	void set_V_time(int beta, int x);
	void display();
	vector<vector<int> > log;
	void check_correctness();

	vector<bool*> ori_U;
	vector<int*> r_U;
	vector<int*> indeg_U;
	vector<int*> neighbor_cnt_U;
	vector<bool*> ori_V;
	vector<int*> r_V;
	vector<int*> indeg_V;
	vector<int*> neighbor_cnt_V;

	void solve_query(int alpha, int beta); vector<int> answer;
};
BD_Index bd;

enum Instruction_Type { INSERT, DELETE, QUERY };
struct Instruction {
	Instruction_Type instruction_type;
	int t1, t2;
	ll time;
};
vector<Instruction> instructions;
vector<double> process_time;
void get_instructions(char* instruction_address) {
	ifstream infile(instruction_address);
	if (!infile.is_open()) {
		check(0, "cannot open file instruction");
		return;
	}

	string line;
	int total_instructions;

	if (getline(infile, line)) {
		stringstream ss(line);
		ss >> total_instructions;
	}

	while (getline(infile, line)) {
		stringstream ss(line);
		char type_char;
		int a, b;
		ll t;

		ss >> type_char >> a >> b >> t;

		Instruction inst;
		if (type_char == 'i') {
			inst.instruction_type = INSERT;
		}
		else if (type_char == 'd') {
			inst.instruction_type = DELETE;
		}
		else if (type_char == 'q') {
			inst.instruction_type = QUERY;
		}
		else {
			check(0, "error instruction type");
		}

		if (type_char == 'i' || type_char == 'd') {
			b += G.U;
			bool found = false;
			for (int j = 0; j < G.undeg[a]; j++) {
				Edge& ne = G.e[G.adj[a][j]];
				if (ne.v == b) {
					a = G.adj[a][j]; found = true;
					break;
				}
			}
			check(found, "not found error");
		}

		inst.t1 = a;
		inst.t2 = b;
		inst.time = t;

		instructions.push_back(inst);
	}
	infile.close();
	process_time.resize(instructions.size());
}

void Graph::read_graph_from_dataset(char* dataset_address) {
	FILE* in = fopen(dataset_address, "r");
	check(in != NULL, "Can not open file dataset_address\n");

	int num1, num2, num3, count;
	char line[256];
	fgets(line, sizeof(line), in);
	count = sscanf(line, "%d %d %d", &num1, &num2, &num3);
	if (count == 2) {
		if (num2 > num1) swap(num1, num2);
		M = num1, U = V = num2;
	}
	else if (count == 3)
		M = num1, U = num2, V = num3;
	else
		check(0, "dataset format error");

	U++, V++; N = U + V;

	e = (Edge*)malloc(M * sizeof(Edge));
	undeg = (int*)malloc(N * sizeof(int)); indeg = (int*)malloc(N * sizeof(int));
	memset(undeg, 0, N * sizeof(int)); memset(indeg, 0, N * sizeof(int));
	adj = (int**)malloc(N * sizeof(int*));
	r = (int*)malloc(N * sizeof(int));

	S.alloc(N); Q.alloc(N); parent.alloc(N); vis.alloc(N); C1.alloc(N); C2.alloc(N); C1_minus_C2.alloc(N); dist.alloc(N); cur.alloc(N); D.alloc(N);

	for (int i = 0; i < M; i++) {
		e[i].u = read_number(in), e[i].v = read_number(in) + U; read_number(in);
		e[i].deleted = false;
		undeg[e[i].u]++, undeg[e[i].v]++;
	}
	d_max = alpha_max = beta_max = 0;
	for (int u = 0; u < U; u++) alpha_max = max(alpha_max, undeg[u]), adj[u] = (int*)malloc(undeg[u] * sizeof(int));
	for (int v = U; v < N; v++) beta_max = max(beta_max, undeg[v]), adj[v] = (int*)malloc(undeg[v] * sizeof(int));
	d_max = max(alpha_max, beta_max) + 2;
	sorted = (int*)malloc(N * sizeof(int)), computed = (bool*)malloc(d_max * sizeof(bool)), position = (int*)malloc(d_max * sizeof(int)), number_of_edges = (int*)malloc(d_max * sizeof(int)), number_of_nodes = (int*)malloc(d_max * sizeof(int));

	memset(undeg, 0, N * sizeof(int));
	for (int i = 0; i < M; i++) {
		int u = e[i].u, v = e[i].v;
		adj[u][undeg[u]++] = i;
		adj[v][undeg[v]++] = i;
	}
}
void Graph::initialize_orientation(int ab, bool the_first_iteration) {
	memset(indeg, 0, N * sizeof(int));
	if (the_first_iteration) {
		for (int u = 0; u < U; u++) {
			if (undeg[u] <= ab) {
				indeg[u] = undeg[u];
				for (int j = 0; j < undeg[u]; j++)
					e[adj[u][j]].to = u;
			}
			else {
				indeg[u] = ab;
				for (int j = 0; j < ab; j++)
					e[adj[u][j]].to = u;
				for (int j = ab; j < undeg[u]; j++)
					e[adj[u][j]].to = e[adj[u][j]].v, indeg[e[adj[u][j]].v]++;
			}
		}
	}
	else {
		for (int v = U; v < N; v++) {
			if (undeg[v] <= ab) {
				indeg[v] = undeg[v];
				for (int j = 0; j < undeg[v]; j++)
					e[adj[v][j]].to = v;
			}
			else {
				indeg[v] = ab;
				for (int j = 0; j < ab; j++)
					e[adj[v][j]].to = v;
				for (int j = ab; j < undeg[v]; j++)
					e[adj[v][j]].to = e[adj[v][j]].u, indeg[e[adj[v][j]].u]++;
			}
		}
	}
}
void Graph::get_core(int core_alpha, int core_beta, Set<int>& C) {
	C.clear();

	// a node --> will_be_deleted_nodes --> deleted_nodes
	Set<int> will_be_deleted_nodes; will_be_deleted_nodes.alloc(N); int pointer = 0;
	Set<int> deleted_nodes; deleted_nodes.alloc(N);
	int* temporary_undeg = (int*)malloc(N * sizeof(int)); memcpy(temporary_undeg, undeg, N * sizeof(int));

	for (int x = 0; x < N; x++)
		if (undeg[x] < (in_U(x) ? core_alpha : core_beta))
			will_be_deleted_nodes.insert(x);

	while (pointer < will_be_deleted_nodes.size) {
		int x = will_be_deleted_nodes.nodes[pointer++];
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			int y = in_U(x) ? ne.v : ne.u;
			if (deleted_nodes.in[y]) continue;
			if (temporary_undeg[y] == (in_U(y) ? core_alpha : core_beta)) {
				will_be_deleted_nodes.insert(y);
			}
			temporary_undeg[x]--;
			temporary_undeg[y]--;
		}
		deleted_nodes.insert(x);
	}

	for (int x = 0; x < N; x++)
		if (!deleted_nodes.in[x])
			C.insert(x);

	free(temporary_undeg);
	return;
}
int max_rank;
void Graph::change_sorted(int position_l, int position_u, Set<int>& D) {
	position_u--;
	while (position_l <= position_u) {
		while (position_l <= position_u && !D.in[sorted[position_l]]) position_l++;
		while (position_l <= position_u && D.in[sorted[position_u]]) position_u--;
		if (position_l <= position_u)
			swap(sorted[position_l], sorted[position_u]);
	}
}
void Graph::Divide(int D_l_rank, int D_u_rank, bool the_first_iteration) {
	// printf("LOG: enter Dividea(%d, %d)\n", D_l_rank, D_u_rank);
	int E_all = number_of_edges[D_l_rank] - number_of_edges[D_u_rank], l = D_l_rank, u = D_u_rank;
	while (u > l) {
		int mid = (u + l + 1) / 2;
		if (the_first_iteration) ReTest(alpha, mid, true);
		else ReTest(mid, beta, false);
		if (number_of_edges[D_l_rank] - number_of_edges[mid] < E_all / 2)
			l = mid;
		else
			u = mid - 1;
	}
	int m = l;
	if (m - D_l_rank >= 2 && number_of_nodes[D_l_rank] != number_of_nodes[m])
		Divide(D_l_rank, m, the_first_iteration);
	m++;
	if (D_u_rank - m >= 2 && number_of_nodes[D_u_rank] != number_of_nodes[m])
		Divide(m, D_u_rank, the_first_iteration);
	return;
}
void Graph::ReTest(int _alpha, int _beta, bool the_first_iteration) {
	alpha = _alpha, beta = _beta;
	// printf("LOG: ReTest(%d, %d)\n", alpha, beta);
	int D_l_rank, D_u_rank;
	if (the_first_iteration) {
		if (computed[beta]) return;
		D_l_rank = beta - 1, D_u_rank = beta + 1;
		while (!computed[D_l_rank]) D_l_rank--;
		while (!computed[D_u_rank]) D_u_rank++;
	}
	else {
		if (computed[alpha]) return;
		D_l_rank = alpha - 1, D_u_rank = alpha + 1;
		while (!computed[D_l_rank]) D_l_rank--;
		while (!computed[D_u_rank]) D_u_rank++;
	}

	C1_minus_C2.clear(); int D_l_position = position[D_l_rank], D_u_position = position[D_u_rank];
	for (int i = D_l_position; i < D_u_position; i++) {
		C1_minus_C2.insert(sorted[i]);
	}

	while (DinicBFS(D_l_rank, D_u_rank)) {
		for (int i = 0; i < C1_minus_C2.size; i++) {
			int x = C1_minus_C2.nodes[i];
			if (in_T(x))
				parent[x] = -2, cur[x] = 0, DinicDFS(x, D_l_rank, D_u_rank);
		}
	}

	int bd = the_first_iteration ? beta : alpha;
	Set<int>& ans = D;
	ans.clear(), Q.clear(), vis.clear();
	for (int i = 0; i < C1_minus_C2.size; i++) {
		int x = C1_minus_C2.nodes[i];
		if (in_T(x))
			Q.push(x), ans.insert(x);
	}
	while (!Q.empty()) {
		int x = Q.pop();
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.to != x) continue;
			int from = in_U(ne.to) ? ne.v : ne.u;
			check(r[from] >= D_l_rank, "r[from] < bd");
			if (r[from] >= D_u_rank) continue;
			if (ans.in[from]) continue;
			Q.push(from), ans.insert(from);
		}
	}
	for (int i = 0; i < ans.size; i++) {
		int x = ans.nodes[i];
		check(r[x] == D_l_rank, "r[x] != D_l_rank");
		r[x] = bd;
	}
	computed[bd] = true;
	number_of_nodes[bd] = ans.size + number_of_nodes[D_u_rank]; position[bd] = D_u_position - ans.size;
	change_sorted(D_l_position, D_u_position, ans);
	int& ans_edges = the_first_iteration ? number_of_edges[beta] : number_of_edges[alpha];
	ans_edges = number_of_edges[D_u_rank];
	for (int i = 0; i < ans.size; i++) ans_edges += indeg[ans.nodes[i]];
	if (ans.size != 0) max_rank = max(max_rank, bd);
	return;
}
void Graph::construct_index() {
	bd.pointer_U.resize(max(alpha_max, beta_max) + 5);
	bd.index_U.resize(max(alpha_max, beta_max) + 5);
	if (algorithm_used == TIME) {
		bd.ori_U.resize(max(alpha_max, beta_max) + 5);
		bd.r_U.resize(max(alpha_max, beta_max) + 5);
		bd.indeg_U.resize(max(alpha_max, beta_max) + 5);
		bd.neighbor_cnt_U.resize(max(alpha_max, beta_max) + 5);
	}

	// printf("----- First iteration -----\n");
	for (alpha = 0; ; alpha++) {
		initialize_orientation(alpha, true);
		max_rank = 0; memset(r, -1, N * sizeof(int));
		for (int i = 0; i < N; i++) sorted[i] = i;
		for (int i = 0; i <= beta_max; i++) computed[i] = false, number_of_edges[i] = 0, number_of_nodes[i] = 0;
		// Get D[0] and D[alpha_max]
		D.clear();
		for (int u = 0; u < U; u++)
			if (undeg[u] > alpha) {
				D.insert(u);
				number_of_edges[0] += undeg[u];
				for (int j = 0; j < undeg[u]; j++) {
					Edge& ne = e[adj[u][j]];
					if (!D.in[ne.v])
						D.insert(ne.v);
				}
			}
		for (int i = 0; i < D.size; i++) {
			int x = D.nodes[i];
			r[x] = 0;
		}
		change_sorted(0, N, D), computed[0] = true, position[0] = N - D.size, number_of_nodes[0] = D.size;
		computed[beta_max] = true, position[beta_max] = N, number_of_edges[beta_max] = 0, number_of_nodes[beta_max] = 0;
		Divide(0, beta_max, true);
		// printf("- %-20s: %d, %d\n", "max non-empty D_{alpha, beta}", alpha, max_rank);
		// analyze_r();
		if (alpha > max_rank) {
			pseudoarboricity = alpha - 1;
			break;
		}

		int nowp = 0; while (r[sorted[nowp]] < alpha) nowp++;
		check(nowp < N, "nowp error");
		int nowr = -1;
		for (int i = nowp; i < N; i++) {
			bd.index_U[alpha].push_back(sorted[i]);
			while (r[sorted[i]] != nowr) nowr++, bd.pointer_U[alpha].push_back(prev(bd.index_U[alpha].end()));
		}

		if (algorithm_used == TIME) {
			bd.ori_U[alpha] = (bool*)malloc(M * sizeof(bool));
			bd.r_U[alpha] = (int*)malloc(N * sizeof(int));
			bd.indeg_U[alpha] = (int*)malloc(N * sizeof(int));
			bd.neighbor_cnt_U[alpha] = (int*)malloc(N * sizeof(int));
			for (int i = 0; i < M; i++) bd.ori_U[alpha][i] = e[i].to == e[i].v;
			for (int x = 0; x < N; x++) bd.r_U[alpha][x] = r[x];
			for (int x = 0; x < N; x++) bd.indeg_U[alpha][x] = indeg[x];
			for (int x = 0; x < N; x++) bd.neighbor_cnt_U[alpha][x] = undeg[x];
		}

		if (alpha == max_rank) {
			pseudoarboricity = alpha;
		}
	}

	bd.pointer_U.resize(pseudoarboricity + 1);
	bd.index_U.resize(pseudoarboricity + 1);
	bd.pointer_V.resize(pseudoarboricity + 1);
	bd.index_V.resize(pseudoarboricity + 1);
	if (algorithm_used == TIME) {
		bd.ori_U.resize(pseudoarboricity + 1);
		bd.r_U.resize(pseudoarboricity + 1);
		bd.indeg_U.resize(pseudoarboricity + 1);
		bd.neighbor_cnt_U.resize(pseudoarboricity + 1);
		bd.ori_V.resize(pseudoarboricity + 1);
		bd.r_V.resize(pseudoarboricity + 1);
		bd.indeg_V.resize(pseudoarboricity + 1);
		bd.neighbor_cnt_V.resize(pseudoarboricity + 1);
	}

	// printf("----- Second iteration -----\n");
	for (beta = 0; ; beta++) {
		initialize_orientation(beta, false);
		max_rank = 0; memset(r, -1, N * sizeof(int));
		for (int i = 0; i < N; i++) sorted[i] = i;
		for (int i = 0; i <= alpha_max; i++) computed[i] = false, number_of_edges[i] = 0, number_of_nodes[i] = 0;
		// Get D[0] and D[alpha_max]
		D.clear();
		for (int v = U; v < N; v++)
			if (undeg[v] > beta) {
				D.insert(v);
				number_of_edges[0] += undeg[v];
				for (int j = 0; j < undeg[v]; j++) {
					Edge& ne = e[adj[v][j]];
					if (!D.in[ne.u])
						D.insert(ne.u);
				}
			}
		for (int i = 0; i < D.size; i++) {
			int x = D.nodes[i];
			r[x] = 0;
		}
		change_sorted(0, N, D), computed[0] = true, position[0] = N - D.size, number_of_nodes[0] = D.size;
		computed[alpha_max] = true, position[alpha_max] = N, number_of_edges[alpha_max] = 0, number_of_nodes[alpha_max] = 0;
		Divide(0, alpha_max, false);
		// printf("- %-20s: %d, %d\n", "max non-empty D_{alpha, beta}", max_rank, beta);
		// analyze_r();
		if (max_rank < beta)
			break;

		int nowp = 0; while (nowp < N && r[sorted[nowp]] <= beta) nowp++;
		int nowr = -1;
		for (int i = nowp; i < N; i++) {
			bd.index_V[beta].push_back(sorted[i]);
			while (r[sorted[i]] != nowr) nowr++, bd.pointer_V[beta].push_back(prev(bd.index_V[beta].end()));
		}

		if (algorithm_used == TIME) {
			bd.ori_V[beta] = (bool*)malloc(M * sizeof(bool));
			bd.r_V[beta] = (int*)malloc(N * sizeof(int));
			bd.indeg_V[beta] = (int*)malloc(N * sizeof(int));
			bd.neighbor_cnt_V[beta] = (int*)malloc(N * sizeof(int));
			for (int i = 0; i < M; i++) bd.ori_V[beta][i] = e[i].to == e[i].v;
			for (int x = 0; x < N; x++) bd.r_V[beta][x] = r[x];
			for (int x = 0; x < N; x++) bd.indeg_V[beta][x] = indeg[x];
			for (int x = 0; x < N; x++) bd.neighbor_cnt_V[beta][x] = undeg[x];
		}

		if (max_rank == beta) // beta = p
			break;
	}

	for (int a = 0; a <= pseudoarboricity; a++)
		for (int b = 0; b < a; b++) bd.pointer_U[a][b] = bd.index_U[a].end();
	for (int b = 0; b <= pseudoarboricity; b++) {
		if (bd.pointer_V[b].size() == 0) bd.pointer_V[b].resize(b + 1);
		for (int a = 0; a <= b; a++) bd.pointer_V[b][a] = bd.index_V[b].end();
	}

	return;
}
bool Graph::DinicBFS(int D_l_rank, int D_u_rank) {
	int dist_t = INF;

	Q.clear(), dist.clear(), parent.clear(), cur.clear();
	for (int i = 0; i < C1_minus_C2.size; i++) {
		int x = C1_minus_C2.nodes[i];
		if (in_T(x))
			dist[x] = 1, Q.push(x);
	}

	bool break_loop = false;
	while (!Q.empty()) {
		int x = Q.pop();
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.to != x) continue;
			int from = in_U(ne.to) ? ne.v : ne.u;
			check(r[from] >= D_l_rank, "r[from] < D_l_rank 2");
			if (r[from] >= D_u_rank) continue;
			if (in_S(from)) {
				dist_t = dist[x] + 2; break_loop = true; break;
			}
			if (dist.in[from]) continue;
			dist[from] = dist[x] + 1;
			Q.push(from);
		}
		if (break_loop) break;
	}
	return dist_t != INF;
}
bool Graph::DinicDFS(int x, int D_l_rank, int D_u_rank) {
	if (in_S(x)) {
		indeg[x]++, indeg[e[parent[x]].to]--, e[parent[x]].to = x;
		return true;
	}
	for (int& j = cur[x]; j < undeg[x]; j++) {
		Edge& ne = e[adj[x][j]];
		if (ne.to != x) continue;
		int from = in_U(ne.to) ? ne.v : ne.u;
		check(r[from] >= D_l_rank, "r[from] < D_l_rank 3");
		if (r[from] >= D_u_rank) continue;
		if ((dist[from] != dist[x] + 1) && !in_S(from)) continue;
		parent[from] = adj[x][j];
		if (DinicDFS(from, D_l_rank, D_u_rank)) {
			if (parent[x] == -2) {
				if (indeg[x] == (in_U(x) ? alpha : beta)) return true;
				continue;
			}
			indeg[x]++, indeg[e[parent[x]].to]--, e[parent[x]].to = x;
			return true;
		}
	}
	return false;
}
bool analyze_r_cmp(int x, int y) {
	return G.r[x] > G.r[y];
}
void Graph::analyze_r() {
	int* sorted_nodes = (int*)malloc(N * sizeof(int));
	for (int i = 0; i < N; i++)	sorted_nodes[i] = i;
	sort(sorted_nodes, sorted_nodes + N, analyze_r_cmp);
	for (int i = 1; i < N; i++) {
		if (r[sorted_nodes[i]] < r[sorted_nodes[i - 1]]) {
			// printf("- %-20s: %d, %d\n", "sorted_nodes[i] size", r[sorted_nodes[i - 1]], i);
		}
	}
	return;
}
void Graph::output_core(Set<int>& C) {
	printf("- %-20s: %d\n", "Core size", C.size);

	if (true) {
		Set<int> output_nodes; output_nodes.alloc(N);
		sort(C.nodes, C.nodes + C.size);
		int i;
		for (i = 0; i < C.size && C.nodes[i] < U; i++)
			output_nodes.insert(C.nodes[i]);
		if (output_nodes.size <= OUTPUT_NODES_LIMIT) {
			printf("- %-20s: ", "Core node in U");
			for (int j = 0; j < output_nodes.size; j++)
				printf("%d ", output_nodes.nodes[j]);
			printf("\n");
		}
		else {
			printf("- %-20s: %d\n", "Core node in U are too many, with size of", output_nodes.size);
		}

		output_nodes.clear();
		for (; i < C.size; i++)
			output_nodes.insert(C.nodes[i]);
		if (output_nodes.size <= OUTPUT_NODES_LIMIT) {
			printf("- %-20s: ", "Core node in V");
			for (int j = 0; j < output_nodes.size; j++)
				printf("%d ", output_nodes.nodes[j]);
			printf("\n");
		}
		else {
			printf("- %-20s: %d\n", "Core node in V are too many, with size of", output_nodes.size);
		}
	}
}
void Graph::output_dense_subgraph() {
	/*
	printf("- %-20s: %d\n", "Dense subgraph size", D.size);

	if (true) {
		Set<int> output_nodes; output_nodes.alloc(N);
		sort(D.nodes, D.nodes + D.size);
		int i;
		for (i = 0; i < D.size && D.nodes[i] < U; i++)
			output_nodes.insert(D.nodes[i]);
		if (output_nodes.size <= OUTPUT_NODES_LIMIT) {
			printf("- %-20s: ", "Dense subgraph node in U");
			for (int j = 0; j < output_nodes.size; j++)
				printf("%d ", output_nodes.nodes[j]);
			printf("\n");
		}
		else {
			printf("- %-20s: %d\n", "Dense subgraph node in U are too many, with size of", output_nodes.size);
		}

		output_nodes.clear();
		for (; i < D.size; i++)
			output_nodes.insert(D.nodes[i]);
		if (output_nodes.size <= OUTPUT_NODES_LIMIT) {
			printf("- %-20s: ", "Dense subgraph node in V");
			for (int j = 0; j < output_nodes.size; j++)
				printf("%d ", output_nodes.nodes[j]);
			printf("\n");
		}
		else {
			printf("- %-20s: %d\n", "Dense subgraph node in V are too many, with size of", output_nodes.size);
		}
	}*/
}
void Graph::display() {
	for (int i = 0; i < M; i++) {
		Edge& ne = e[i];
		int from = in_U(ne.to) ? ne.v : ne.u;
		printf("%d %d\n", from, ne.to);
	}
}
void Graph::check_correctness() {
	for (int i = 0; i < D.size; i++) {
		printf("%d, ", D.nodes[i]);
	}
	printf("\n");
	return;
}
void Graph::Dynamic_ReTest(int _alpha, int _beta) {
	alpha = _alpha, beta = _beta;

	C1_minus_C2.clear();
	for (int i = 0; i < C1.size; i++) {
		int x = C1.nodes[i];
		if (!C2.in[x])
			C1_minus_C2.insert(x);
	}

	for (int i = 0; i < C1_minus_C2.size; i++) {
		indeg[C1_minus_C2.nodes[i]] = 0;
	}
	for (int i = 0; i < C1_minus_C2.size; i++) {
		int x = C1_minus_C2.nodes[i];
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.deleted) continue;
			int y = in_U(x) ? ne.v : ne.u;
			if (C2.in[y]) {
				ne.to = x, indeg[ne.to]++;
				continue;
			}
			else if (!C1.in[y]) {
				ne.to = y;
				continue;
			}
			else {
				if (y == ne.v) ne.to = y, indeg[ne.to]++;
				continue;
			}
		}
	}

	while (Dynamic_DinicBFS()) {
		for (int i = 0; i < C1_minus_C2.size; i++) {
			int x = C1_minus_C2.nodes[i];
			if (in_T(x))
				parent[x] = -2, cur[x] = 0, Dynamic_DinicDFS(x);
		}
	}

	D.clear(), Q.clear(), vis.clear();
	for (int i = 0; i < C1_minus_C2.size; i++) {
		int x = C1_minus_C2.nodes[i];
		if (in_T(x))
			Q.push(x), D.insert(x);
	}
	while (!Q.empty()) {
		int x = Q.pop();
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.deleted) continue;
			if (ne.to != x) continue;
			int from = in_U(ne.to) ? ne.v : ne.u;
			check(C1.in[from], "C1_minus_C2 error");
			if (C2.in[from]) continue;
			if (D.in[from]) continue;
			Q.push(from), D.insert(from);
		}
	}
	for (int i = 0; i < C2.size; i++) {
		int x = C2.nodes[i];
		check(!D.in[x], "C2 error");
		D.insert(x);
	}
	return;
}
bool Graph::Dynamic_DinicBFS() {
	int dist_t = INF;

	Q.clear(), dist.clear(), parent.clear(), cur.clear();
	for (int i = 0; i < C1_minus_C2.size; i++) {
		int x = C1_minus_C2.nodes[i];
		if (in_T(x))
			dist[x] = 1, Q.push(x);
	}

	bool break_loop = false;
	while (!Q.empty()) {
		int x = Q.pop();
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.deleted) continue;
			if (ne.to != x) continue;
			int from = in_U(ne.to) ? ne.v : ne.u;
			check(C1.in[from], "C1_minus_C2 error");
			if (C2.in[from]) continue;
			if (in_S(from)) {
				dist_t = dist[x] + 2; break_loop = true; break;
			}
			if (dist.in[from]) continue;
			dist[from] = dist[x] + 1;
			Q.push(from);
		}
		if (break_loop) break;
	}
	return dist_t != INF;
}
bool Graph::Dynamic_DinicDFS(int x) {
	if (in_S(x)) {
		indeg[x]++, indeg[e[parent[x]].to]--, e[parent[x]].to = x;
		return true;
	}
	for (int& j = cur[x]; j < undeg[x]; j++) {
		Edge& ne = e[adj[x][j]];
		if (ne.deleted) continue;
		if (ne.to != x) continue;
		int from = in_U(ne.to) ? ne.v : ne.u;
		check(C1.in[from], "C1_minus_C2 error");
		if (C2.in[from]) continue;
		if ((dist[from] != dist[x] + 1) && !in_S(from)) continue;
		parent[from] = adj[x][j];
		if (Dynamic_DinicDFS(from)) {
			if (parent[x] == -2) {
				if (indeg[x] == (in_U(x) ? alpha : beta)) return true;
				continue;
			}
			indeg[x]++, indeg[e[parent[x]].to]--, e[parent[x]].to = x;
			return true;
		}
	}
	return false;
}
void Graph::delete_edge_space(int edge_id) {
	int u = e[edge_id].u, v = e[edge_id].v;
	e[edge_id].deleted = true;
	unordered_map<int, int> neighbor;
	for (int j = 0; j < undeg[u]; j++) {
		Edge& ne = e[adj[u][j]];
		if (ne.deleted && adj[u][j] != edge_id) continue;
		neighbor[ne.v] = -1;
	}
	C1.clear(), C2.clear();
	for (int x = 0; x < N; x++) C1.insert(x);
	for (int a = 0; a <= pseudoarboricity; a++) {
		for (auto p : neighbor) neighbor[p.first] = a - 1;
		int len = bd.pointer_U[a].size(), u_rank = a - 1, v_rank = a - 1;
		for (int b = a; b < len; b++) {
			list<int>::iterator it_start = bd.pointer_U[a][b];
			list<int>::iterator it_end = b == len - 1 ? bd.index_U[a].end() : bd.pointer_U[a][b + 1];
			while (it_start != it_end) {
				if (neighbor.count(*it_start)) neighbor[*it_start] = b;
				if (*it_start == u) u_rank = b;
				if (*it_start == v) v_rank = b;
				it_start++;
			}
		}
		if (u_rank > v_rank) {
			if (v_rank != a - 1) {
				Dynamic_ReTest(a, v_rank);
				unordered_set<int> update_nodes;
				list<int>::iterator it_start = bd.pointer_U[a][v_rank];
				list<int>::iterator it_end = v_rank == len - 1 ? bd.index_U[a].end() : bd.pointer_U[a][v_rank + 1];
				while (it_start != it_end) {
					if (!D.in[*it_start]) update_nodes.insert(*it_start);
					it_start++;
				}
				bd.dec_U(a, v_rank, update_nodes);
			}
		}
		else if (u_rank <= v_rank) {
			vector<int> sorted_neighbor;
			for (auto p : neighbor) sorted_neighbor.push_back(p.second);
			sort(sorted_neighbor.begin(), sorted_neighbor.end(), greater<int>());
			int nowr = sorted_neighbor.size() > a ? sorted_neighbor[a] : a - 1;
			if (nowr != a - 1) {
				Dynamic_ReTest(a, nowr);
				unordered_set<int> update_nodes;
				list<int>::iterator it_start = bd.pointer_U[a][nowr];
				list<int>::iterator it_end = nowr == len - 1 ? bd.index_U[a].end() : bd.pointer_U[a][nowr + 1];
				while (it_start != it_end) {
					if (!D.in[*it_start]) update_nodes.insert(*it_start);
					it_start++;
				}
				bd.dec_U(a, nowr, update_nodes);
			}
			bd.set_U(a, u);
		}
	}

	neighbor.clear();
	for (int j = 0; j < undeg[v]; j++) {
		Edge& ne = e[adj[v][j]];
		if (ne.deleted && adj[v][j] != edge_id) continue;
		neighbor[ne.u] = -1;
	}
	for (int b = 0; b <= pseudoarboricity; b++) {
		for (auto p : neighbor) neighbor[p.first] = b;
		int len = bd.pointer_V[b].size(), u_rank = b, v_rank = b;
		for (int a = b + 1; a < len; a++) {
			list<int>::iterator it_start = bd.pointer_V[b][a];
			list<int>::iterator it_end = a == len - 1 ? bd.index_V[b].end() : bd.pointer_V[b][a + 1];
			while (it_start != it_end) {
				if (neighbor.count(*it_start)) neighbor[*it_start] = a;
				if (*it_start == v) v_rank = a;
				if (*it_start == u) u_rank = a;
				it_start++;
			}
		}
		if (v_rank > u_rank) {
			if (u_rank != b) {
				Dynamic_ReTest(u_rank, b);
				unordered_set<int> update_nodes;
				list<int>::iterator it_start = bd.pointer_V[b][u_rank];
				list<int>::iterator it_end = u_rank == len - 1 ? bd.index_V[b].end() : bd.pointer_V[b][u_rank + 1];
				while (it_start != it_end) {
					if (!D.in[*it_start]) update_nodes.insert(*it_start);
					it_start++;
				}
				bd.dec_V(b, u_rank, update_nodes);
			}
		}
		else if (v_rank <= u_rank) {
			vector<int> sorted_neighbor;
			for (auto p : neighbor) sorted_neighbor.push_back(p.second);
			sort(sorted_neighbor.begin(), sorted_neighbor.end(), greater<int>());
			int nowr = sorted_neighbor.size() > b ? sorted_neighbor[b] : b;
			if (nowr != b) {
				Dynamic_ReTest(nowr, b);
				unordered_set<int> update_nodes;
				list<int>::iterator it_start = bd.pointer_V[b][nowr];
				list<int>::iterator it_end = nowr == len - 1 ? bd.index_V[b].end() : bd.pointer_V[b][nowr + 1];
				while (it_start != it_end) {
					if (!D.in[*it_start]) update_nodes.insert(*it_start);
					it_start++;
				}
				bd.dec_V(b, nowr, update_nodes);
			}
			bd.set_V(b, v);
		}
	}
}
void Graph::insert_edge_space(int edge_id) {
	int u = e[edge_id].u, v = e[edge_id].v;
	e[edge_id].deleted = false;
	unordered_map<int, int> neighbor;
	for (int j = 0; j < undeg[u]; j++) {
		Edge& ne = e[adj[u][j]];
		if (ne.deleted) continue;
		neighbor[ne.v] = -1;
	}
	C1.clear(), C2.clear();
	for (int x = 0; x < N; x++) C1.insert(x);
	for (int a = 0; a <= pseudoarboricity; a++) {
		for (auto p : neighbor) neighbor[p.first] = a - 1;
		int len = bd.pointer_U[a].size(), u_rank = a - 1, v_rank = a - 1;
		for (int b = a; b < len; b++) {
			list<int>::iterator it_start = bd.pointer_U[a][b];
			list<int>::iterator it_end = b == len - 1 ? bd.index_U[a].end() : bd.pointer_U[a][b + 1];
			while (it_start != it_end) {
				if (neighbor.count(*it_start)) neighbor[*it_start] = b;
				if (*it_start == u) u_rank = b;
				if (*it_start == v) v_rank = b;
				it_start++;
			}
		}
		if (u_rank > v_rank) {
			Dynamic_ReTest(a, v_rank + 1);
			unordered_set<int> update_nodes;
			if (v_rank > alpha - 1) {
				list<int>::iterator it_start = bd.pointer_U[a][v_rank];
				list<int>::iterator it_end = v_rank == len - 1 ? bd.index_U[a].end() : bd.pointer_U[a][v_rank + 1];
				while (it_start != it_end) {
					if (D.in[*it_start]) update_nodes.insert(*it_start);
					it_start++;
				}
			}
			else {
				for (int i = 0; i < D.size; i++) update_nodes.insert(D.nodes[i]);
				list<int>::iterator it_start = bd.index_U[a].begin();
				while (it_start != bd.index_U[a].end()) update_nodes.erase(*it_start), it_start++;
			}
			bd.inc_U(a, v_rank, update_nodes);
		}
		else if (u_rank <= v_rank) {
			vector<int> sorted_neighbor;
			for (auto p : neighbor) sorted_neighbor.push_back(p.second);
			sort(sorted_neighbor.begin(), sorted_neighbor.end(), greater<int>());
			int nowr = sorted_neighbor.size() > a ? sorted_neighbor[a] : a - 1;
			Dynamic_ReTest(a, nowr + 1);
			unordered_set<int> update_nodes;
			if (nowr > alpha - 1) {
				list<int>::iterator it_start = bd.pointer_U[a][nowr];
				list<int>::iterator it_end = nowr == len - 1 ? bd.index_U[a].end() : bd.pointer_U[a][nowr + 1];
				while (it_start != it_end) {
					if (D.in[*it_start]) update_nodes.insert(*it_start);
					it_start++;
				}
			}
			else {
				for (int i = 0; i < D.size; i++) update_nodes.insert(D.nodes[i]);
				list<int>::iterator it_start = bd.index_U[a].begin();
				while (it_start != bd.index_U[a].end()) update_nodes.erase(*it_start), it_start++;
			}
			bd.inc_U(a, nowr, update_nodes);
			bd.set_U(a, u);
		}
	}

	neighbor.clear();
	for (int j = 0; j < undeg[v]; j++) {
		Edge& ne = e[adj[v][j]];
		if (ne.deleted) continue;
		neighbor[ne.u] = -1;
	}
	for (int b = 0; b <= pseudoarboricity; b++) {
		for (auto p : neighbor) neighbor[p.first] = b;
		int len = bd.pointer_V[b].size(), u_rank = b, v_rank = b;
		for (int a = b + 1; a < len; a++) {
			list<int>::iterator it_start = bd.pointer_V[b][a];
			list<int>::iterator it_end = a == len - 1 ? bd.index_V[b].end() : bd.pointer_V[b][a + 1];
			while (it_start != it_end) {
				if (neighbor.count(*it_start)) neighbor[*it_start] = a;
				if (*it_start == u) u_rank = a;
				if (*it_start == v) v_rank = a;
				it_start++;
			}
		}
		if (v_rank > u_rank) {
			Dynamic_ReTest(u_rank + 1, b);
			unordered_set<int> update_nodes;
			if (u_rank > beta) {
				list<int>::iterator it_start = bd.pointer_V[b][u_rank];
				list<int>::iterator it_end = u_rank == len - 1 ? bd.index_V[b].end() : bd.pointer_V[b][u_rank + 1];
				while (it_start != it_end) {
					if (D.in[*it_start]) update_nodes.insert(*it_start);
					it_start++;
				}
			}
			else {
				for (int i = 0; i < D.size; i++) update_nodes.insert(D.nodes[i]);
				list<int>::iterator it_start = bd.index_V[b].begin();
				while (it_start != bd.index_V[b].end()) update_nodes.erase(*it_start), it_start++;
			}
			bd.inc_V(b, u_rank, update_nodes);
		}
		else if (v_rank <= u_rank) {
			vector<int> sorted_neighbor;
			for (auto p : neighbor) sorted_neighbor.push_back(p.second);
			sort(sorted_neighbor.begin(), sorted_neighbor.end(), greater<int>());
			int nowr = sorted_neighbor.size() > b ? sorted_neighbor[b] : b;
			Dynamic_ReTest(nowr + 1, b);
			unordered_set<int> update_nodes;
			if (nowr > beta) {
				list<int>::iterator it_start = bd.pointer_V[b][nowr];
				list<int>::iterator it_end = nowr == len - 1 ? bd.index_V[b].end() : bd.pointer_V[b][nowr + 1];
				while (it_start != it_end) {
					if (D.in[*it_start]) update_nodes.insert(*it_start);
					it_start++;
				}
			}
			else {
				for (int i = 0; i < D.size; i++) update_nodes.insert(D.nodes[i]);
				list<int>::iterator it_start = bd.index_V[b].begin();
				while (it_start != bd.index_V[b].end()) update_nodes.erase(*it_start), it_start++;
			}
			bd.inc_V(b, nowr, update_nodes);
			bd.set_V(b, v);
		}
	}
}
void Graph::delete_edge_time(int edge_id) {
	int u = e[edge_id].u, v = e[edge_id].v;
	for (int a = 0; a <= pseudoarboricity; a++) {
		e[edge_id].deleted = false;
		if (!bd.ori_U[a][edge_id]) {
			if (bd.neighbor_cnt_U[a][u] <= a) {
				e[edge_id].deleted = true, bd.indeg_U[a][u]--;
			}
			else {
				int r0 = bd.r_U[a][u], w = -1; bool have_finished = false;
				Q.clear(), parent.clear();
				Q.push(u), parent[u] = -1;
				while (!Q.empty()) {
					int x = Q.pop();
					for (int j = 0; j < undeg[x]; j++) {
						int ne = adj[x][j];
						if (e[ne].deleted) continue;
						if (in_V(x) == bd.ori_U[a][ne]) continue;
						int to = in_U(x) ? e[ne].v : e[ne].u;
						if (bd.r_U[a][to] != r0) continue;
						if (parent.in[to]) continue;
						if (in_V(to)) {
							if (bd.indeg_U[a][to] == r0 + 1) {
								parent[to] = ne, w = to, have_finished = true;
								break;
							}
						}
						parent[to] = ne, Q.push(to);
					}
					if (have_finished) break;
				}
				check(w != -1, "w error");
				while (w != u) {
					int from, to;
					if (bd.ori_U[a][parent[w]]) from = e[parent[w]].u, to = e[parent[w]].v;
					else from = e[parent[w]].v, to = e[parent[w]].u;
					bd.indeg_U[a][to]--, bd.indeg_U[a][from]++;
					bd.ori_U[a][parent[w]] = !bd.ori_U[a][parent[w]];
					w = from;
				}
				e[edge_id].deleted = true, bd.indeg_U[a][u]--;
				check(!bd.ori_U[a][edge_id], "bd.ori_U[a][edge_id] error");

				Q.clear(), D.clear(), vis.clear();
				for (int u = 0; u < U; u++)
					if (bd.r_U[a][u] == r0)
						vis.insert(u);
				for (int v = U; v < N; v++) {
					if (bd.r_U[a][v] == r0) {
						vis.insert(v);
						if (bd.indeg_U[a][v] == r0 + 1)
							Q.push(v), D.insert(v);
					}
				}
				while (!Q.empty()) {
					int x = Q.pop();
					for (int j = 0; j < undeg[x]; j++) {
						int ne = adj[x][j];
						if (e[ne].deleted) continue;
						if (in_U(x) == bd.ori_U[a][ne]) continue;
						int from = in_U(x) ? e[ne].v : e[ne].u;
						if (bd.r_U[a][from] != r0) continue;
						if (D.in[from]) continue;
						Q.push(from), D.insert(from);
					}
				}
				unordered_set<int> update_nodes;
				for (int i = 0; i < vis.size; i++)
					if (!D.in[vis.nodes[i]]) update_nodes.insert(vis.nodes[i]), bd.r_U[a][vis.nodes[i]]--;
				if (r0 >= a)
					bd.dec_U(a, r0, update_nodes);
				bd.set_U_time(a, u);
			}
		}
		else {
			int r0 = bd.r_U[a][v], w = -1; bool have_finished = false;
			if (bd.indeg_U[a][v] == r0 + 1) {
				w = v;
				goto BFS_end;
			}
			Q.clear(), parent.clear();
			Q.push(v), parent[v] = -1;
			while (!Q.empty()) {
				int x = Q.pop();
				for (int j = 0; j < undeg[x]; j++) {
					int ne = adj[x][j];
					if (e[ne].deleted) continue;
					if (in_V(x) == bd.ori_U[a][ne]) continue;
					int to = in_U(x) ? e[ne].v : e[ne].u;
					if (bd.r_U[a][to] != r0) continue;
					if (parent.in[to]) continue;
					if (in_V(to)) {
						if (bd.indeg_U[a][to] == r0 + 1) {
							parent[to] = ne, w = to, have_finished = true;
							break;
						}
					}
					parent[to] = ne, Q.push(to);
				}
				if (have_finished) break;
			}
		BFS_end:
			check(w != -1, "w error");
			while (w != v) {
				int from, to;
				if (bd.ori_U[a][parent[w]]) from = e[parent[w]].u, to = e[parent[w]].v;
				else from = e[parent[w]].v, to = e[parent[w]].u;
				bd.indeg_U[a][to]--, bd.indeg_U[a][from]++;
				bd.ori_U[a][parent[w]] = !bd.ori_U[a][parent[w]];
				w = from;
			}
			e[edge_id].deleted = true, bd.indeg_U[a][v]--;
			check(bd.ori_U[a][edge_id], "bd.ori_U[a][edge_id] error 2");

			Q.clear(), D.clear(), vis.clear();
			for (int u = 0; u < U; u++)
				if (bd.r_U[a][u] == r0)
					vis.insert(u);
			for (int v = U; v < N; v++) {
				if (bd.r_U[a][v] == r0) {
					vis.insert(v);
					if (bd.indeg_U[a][v] == r0 + 1)
						Q.push(v), D.insert(v);
				}
			}
			while (!Q.empty()) {
				int x = Q.pop();
				for (int j = 0; j < undeg[x]; j++) {
					int ne = adj[x][j];
					if (e[ne].deleted) continue;
					if (in_U(x) == bd.ori_U[a][ne]) continue;
					int from = in_U(x) ? e[ne].v : e[ne].u;
					if (bd.r_U[a][from] != r0) continue;
					if (D.in[from]) continue;
					Q.push(from), D.insert(from);
				}
			}
			unordered_set<int> update_nodes;
			for (int i = 0; i < vis.size; i++)
				if (!D.in[vis.nodes[i]]) update_nodes.insert(vis.nodes[i]), bd.r_U[a][vis.nodes[i]]--;
			if (r0 >= a)
				bd.dec_U(a, r0, update_nodes);
			bd.set_U_time(a, u);
		}
	}
	e[edge_id].deleted = true;
	for (int a = 0; a <= pseudoarboricity; a++)
		bd.neighbor_cnt_U[a][e[edge_id].u]--, bd.neighbor_cnt_U[a][e[edge_id].v]--;

	for (int b = 0; b <= pseudoarboricity; b++) {
		e[edge_id].deleted = false;
		if (bd.ori_V[b][edge_id]) {
			if (bd.neighbor_cnt_V[b][v] <= b) {
				e[edge_id].deleted = true, bd.indeg_V[b][v]--;
			}
			else {
				int r0 = bd.r_V[b][v], w = -1; bool have_finished = false;
				Q.clear(), parent.clear();
				Q.push(v), parent[v] = -1;
				while (!Q.empty()) {
					int x = Q.pop();
					for (int j = 0; j < undeg[x]; j++) {
						int ne = adj[x][j];
						if (e[ne].deleted) continue;
						if (in_V(x) == bd.ori_V[b][ne]) continue;
						int to = in_U(x) ? e[ne].v : e[ne].u;
						if (bd.r_V[b][to] != r0) continue;
						if (parent.in[to]) continue;
						if (in_U(to)) {
							if (bd.indeg_V[b][to] == r0 + 1) {
								parent[to] = ne, w = to, have_finished = true;
								break;
							}
						}
						parent[to] = ne, Q.push(to);
					}
					if (have_finished) break;
				}
				check(w != -1, "w error");
				while (w != v) {
					int from, to;
					if (bd.ori_V[b][parent[w]]) from = e[parent[w]].u, to = e[parent[w]].v;
					else from = e[parent[w]].v, to = e[parent[w]].u;
					bd.indeg_V[b][to]--, bd.indeg_V[b][from]++;
					bd.ori_V[b][parent[w]] = !bd.ori_V[b][parent[w]];
					w = from;
				}
				e[edge_id].deleted = true, bd.indeg_V[b][v]--;
				check(bd.ori_V[b][edge_id], "bd.ori_V[b][edge_id] error");

				Q.clear(), D.clear(), vis.clear();
				for (int v = U; v < N; v++)
					if (bd.r_V[b][v] == r0)
						vis.insert(v);
				for (int u = 0; u < U; u++) {
					if (bd.r_V[b][u] == r0) {
						vis.insert(u);
						if (bd.indeg_V[b][u] == r0 + 1)
							Q.push(u), D.insert(u);
					}
				}
				while (!Q.empty()) {
					int x = Q.pop();
					for (int j = 0; j < undeg[x]; j++) {
						int ne = adj[x][j];
						if (e[ne].deleted) continue;
						if (in_U(x) == bd.ori_V[b][ne]) continue;
						int from = in_U(x) ? e[ne].v : e[ne].u;
						if (bd.r_V[b][from] != r0) continue;
						if (D.in[from]) continue;
						Q.push(from), D.insert(from);
					}
				}
				unordered_set<int> update_nodes;
				for (int i = 0; i < vis.size; i++)
					if (!D.in[vis.nodes[i]]) update_nodes.insert(vis.nodes[i]), bd.r_V[b][vis.nodes[i]]--;
				if (r0 >= b + 1)
					bd.dec_V(b, r0, update_nodes);
				bd.set_V_time(b, v);
			}
		}
		else {
			int r0 = bd.r_V[b][u], w = -1; bool have_finished = false;
			if (bd.indeg_V[b][u] == r0 + 1) {
				w = u;
				goto BFS_end_V;
			}
			Q.clear(), parent.clear();
			Q.push(u), parent[u] = -1;
			while (!Q.empty()) {
				int x = Q.pop();
				for (int j = 0; j < undeg[x]; j++) {
					int ne = adj[x][j];
					if (e[ne].deleted) continue;
					if (in_V(x) == bd.ori_V[b][ne]) continue;
					int to = in_U(x) ? e[ne].v : e[ne].u;
					if (bd.r_V[b][to] != r0) continue;
					if (parent.in[to]) continue;
					if (in_U(to)) {
						if (bd.indeg_V[b][to] == r0 + 1) {
							parent[to] = ne, w = to, have_finished = true;
							break;
						}
					}
					parent[to] = ne, Q.push(to);
				}
				if (have_finished) break;
			}
		BFS_end_V:
			check(w != -1, "w error");
			while (w != u) {
				int from, to;
				if (bd.ori_V[b][parent[w]]) from = e[parent[w]].u, to = e[parent[w]].v;
				else from = e[parent[w]].v, to = e[parent[w]].u;
				bd.indeg_V[b][to]--, bd.indeg_V[b][from]++;
				bd.ori_V[b][parent[w]] = !bd.ori_V[b][parent[w]];
				w = from;
			}
			e[edge_id].deleted = true, bd.indeg_V[b][u]--;
			check(!bd.ori_V[b][edge_id], "bd.ori_V[b][edge_id] error 2");

			Q.clear(), D.clear(), vis.clear();
			for (int v = U; v < N; v++)
				if (bd.r_V[b][v] == r0)
					vis.insert(v);
			for (int u = 0; u < U; u++) {
				if (bd.r_V[b][u] == r0) {
					vis.insert(u);
					if (bd.indeg_V[b][u] == r0 + 1)
						Q.push(u), D.insert(u);
				}
			}
			while (!Q.empty()) {
				int x = Q.pop();
				for (int j = 0; j < undeg[x]; j++) {
					int ne = adj[x][j];
					if (e[ne].deleted) continue;
					if (in_U(x) == bd.ori_V[b][ne]) continue;
					int from = in_U(x) ? e[ne].v : e[ne].u;
					if (bd.r_V[b][from] != r0) continue;
					if (D.in[from]) continue;
					Q.push(from), D.insert(from);
				}
			}
			unordered_set<int> update_nodes;
			for (int i = 0; i < vis.size; i++)
				if (!D.in[vis.nodes[i]]) update_nodes.insert(vis.nodes[i]), bd.r_V[b][vis.nodes[i]]--;
			if (r0 >= b + 1)
				bd.dec_V(b, r0, update_nodes);
			bd.set_V_time(b, v);
		}
	}
	e[edge_id].deleted = true;
	for (int b = 0; b <= pseudoarboricity; b++)
		bd.neighbor_cnt_V[b][e[edge_id].u]--, bd.neighbor_cnt_V[b][e[edge_id].v]--;
}
void Graph::insert_edge_time(int edge_id) {
	int u = e[edge_id].u, v = e[edge_id].v;
	for (int a = 0; a <= pseudoarboricity; a++) {
		e[edge_id].deleted = true;
		int r0_u = bd.r_U[a][u], r0_v = bd.r_U[a][v];
		bd.ori_U[a][edge_id] = false;
		bd.indeg_U[a][u]++;
		e[edge_id].deleted = false;
		if (bd.neighbor_cnt_U[a][u] < a) continue;

		Q.clear(), parent.clear();
		int r0 = INF;
		for (int j = 0; j < undeg[u]; j++) {
			int ne = adj[u][j];
			if (e[ne].deleted) continue;
			if (bd.ori_U[a][ne]) continue;
			if (bd.r_U[a][e[ne].v] < r0) {
				r0 = bd.r_U[a][e[ne].v];
				parent.clear(); parent[e[ne].v] = ne;
			}
			else if (bd.r_U[a][e[ne].v] == r0) {
				parent[e[ne].v] = ne;
			}
		}
		if (r0 == INF) r0 = -1;
		int w = -1;
		bool have_finished = true;
		for (int i = 0; i < parent.size; i++) {
			int v = parent.nodes[i];
			if (bd.indeg_U[a][v] == r0) {
				w = v;
				goto BFS_end;
			}
		}
		parent[u] = -1;
		for (int i = 0; i < parent.size; i++) Q.push(parent.nodes[i]);
		while (!Q.empty()) {
			int x = Q.pop();
			for (int j = 0; j < undeg[x]; j++) {
				int ne = adj[x][j];
				if (e[ne].deleted) continue;
				if (in_U(x) == bd.ori_U[a][ne]) continue;
				int from = in_U(x) ? e[ne].v : e[ne].u;
				if (bd.r_U[a][from] != r0) continue;
				if (parent.in[from]) continue;
				if (in_V(from) && bd.indeg_U[a][from] == r0) {
					parent[from] = ne, w = from, have_finished = false;
					break;
				}
				parent[from] = ne, Q.push(from);
			}
		}
	BFS_end:
		if (w == -1) {
			w = parent.nodes[0];
			bd.indeg_U[a][u]--;
			bd.indeg_U[a][w]++;
			bd.ori_U[a][parent[w]] = true;
			Q.clear(), Q.push(w);
			unordered_set<int> update_nodes;
			bd.r_U[a][w] = r0 + 1, update_nodes.insert(w);
			while (!Q.empty()) {
				int x = Q.pop();
				for (int j = 0; j < undeg[x]; j++) {
					int ne = adj[x][j];
					if (e[ne].deleted) continue;
					if (in_U(x) == bd.ori_U[a][ne]) continue;
					int from = in_U(x) ? e[ne].v : e[ne].u;
					check(bd.r_U[a][from] >= r0 || u == from, "bd.r_U[a][from] >= r0 error");
					if (bd.r_U[a][from] != r0 && u != from) continue;
					if (from != u)
						bd.r_U[a][from] = r0 + 1, update_nodes.insert(from);
					Q.push(from);
				}
			}
			if (r0 >= a - 1)
				bd.inc_U(a, r0, update_nodes);
		}
		else {
			while (w != u) {
				int from, to;
				if (bd.ori_U[a][parent[w]]) from = e[parent[w]].u, to = e[parent[w]].v;
				else from = e[parent[w]].v, to = e[parent[w]].u;
				bd.indeg_U[a][to]--, bd.indeg_U[a][from]++;
				bd.ori_U[a][parent[w]] = !bd.ori_U[a][parent[w]];
				w = to;
			}
			check(bd.indeg_U[a][u] == a, "bd.indeg_U[a][u] == a error");
		}
		if (r0_u <= r0_v) {
			bd.set_U_time(a, u);
		}
	}
	e[edge_id].deleted = false;
	for (int a = 0; a <= pseudoarboricity; a++)
		bd.neighbor_cnt_U[a][e[edge_id].u]++, bd.neighbor_cnt_U[a][e[edge_id].v]++;

	for (int b = 0; b <= pseudoarboricity; b++) {
		e[edge_id].deleted = true;
		int r0_v = bd.r_V[b][v], r0_u = bd.r_V[b][u];
		bd.ori_V[b][edge_id] = true;
		bd.indeg_V[b][v]++;
		e[edge_id].deleted = false;
		if (bd.neighbor_cnt_V[b][v] < b) continue;

		Q.clear(), parent.clear();
		int r0 = INF;
		for (int j = 0; j < undeg[v]; j++) {
			int ne = adj[v][j];
			if (e[ne].deleted) continue;
			if (!bd.ori_V[b][ne]) continue;
			if (bd.r_V[b][e[ne].u] < r0) {
				r0 = bd.r_V[b][e[ne].u];
				parent.clear(); parent[e[ne].u] = ne;
			}
			else if (bd.r_V[b][e[ne].u] == r0) {
				parent[e[ne].u] = ne;
			}
		}
		if (r0 == INF) r0 = -1;
		int w = -1;
		bool have_finished = true;
		for (int i = 0; i < parent.size; i++) {
			int u = parent.nodes[i];
			if (bd.indeg_V[b][u] == r0) {
				w = u;
				goto BFS_end_V;
			}
		}
		parent[v] = -1;
		for (int i = 0; i < parent.size; i++) Q.push(parent.nodes[i]);
		while (!Q.empty()) {
			int x = Q.pop();
			for (int j = 0; j < undeg[x]; j++) {
				int ne = adj[x][j];
				if (e[ne].deleted) continue;
				if (in_U(x) == bd.ori_V[b][ne]) continue;
				int from = in_U(x) ? e[ne].v : e[ne].u;
				if (bd.r_V[b][from] != r0) continue;
				if (parent.in[from]) continue;
				if (in_U(from) && bd.indeg_V[b][from] == r0) {
					parent[from] = ne, w = from, have_finished = false;
					break;
				}
				parent[from] = ne, Q.push(from);
			}
		}
	BFS_end_V:
		if (w == -1) {
			w = parent.nodes[0];
			bd.indeg_V[b][v]--;
			bd.indeg_V[b][w]++;
			bd.ori_V[b][parent[w]] = false;
			Q.clear(), Q.push(w);
			unordered_set<int> update_nodes;
			bd.r_V[b][w] = r0 + 1, update_nodes.insert(w);
			while (!Q.empty()) {
				int x = Q.pop();
				for (int j = 0; j < undeg[x]; j++) {
					int ne = adj[x][j];
					if (e[ne].deleted) continue;
					if (in_U(x) == bd.ori_V[b][ne]) continue;
					int from = in_U(x) ? e[ne].v : e[ne].u;
					check(bd.r_V[b][from] >= r0 || v == from, "bd.r_V[b][from] >= r0 error");
					if (bd.r_V[b][from] != r0 && v != from) continue;
					if (from != v)
						bd.r_V[b][from] = r0 + 1, update_nodes.insert(from);
					Q.push(from);
				}
			}
			if (r0 >= b)
				bd.inc_V(b, r0, update_nodes);
		}
		else {
			while (w != v) {
				int from, to;
				if (bd.ori_V[b][parent[w]]) from = e[parent[w]].u, to = e[parent[w]].v;
				else from = e[parent[w]].v, to = e[parent[w]].u;
				bd.indeg_V[b][to]--, bd.indeg_V[b][from]++;
				bd.ori_V[b][parent[w]] = !bd.ori_V[b][parent[w]];
				w = to;
			}
			check(bd.indeg_V[b][v] == b, "bd.indeg_U[b][v] == b error");
		}
		if (r0_v <= r0_u) {
			bd.set_V_time(b, v);
		}
	}
	e[edge_id].deleted = false;
	for (int b = 0; b <= pseudoarboricity; b++)
		bd.neighbor_cnt_V[b][e[edge_id].u]++, bd.neighbor_cnt_V[b][e[edge_id].v]++;
}

void BD_Index::inc_U(int alpha, int beta, unordered_set<int>& update_nodes) {
	if (update_nodes.empty()) return;
	if (beta == alpha - 1) {
		for (auto x : update_nodes)
			index_U[alpha].push_front(x);
		pointer_U[alpha][beta + 1] = index_U[alpha].begin();
		return;
	}
	list<int>::iterator it_start = pointer_U[alpha][beta];
	check(pointer_U[alpha].size() > beta + 1, "pointer_U[alpha].size() error");
	list<int>::iterator it_end = pointer_U[alpha][beta + 1];
	int start_beta = beta; while (start_beta >= 0 && pointer_U[alpha][start_beta] == pointer_U[alpha][beta]) start_beta--;

	unordered_set<int> no_update_nodes;
	while (it_start != it_end) {
		if (!update_nodes.count(*it_start))
			no_update_nodes.insert(*it_start);
		it_start = index_U[alpha].erase(it_start);
	}
	bool beta_pointer = false, beta_plus_pointer = false;
	for (auto x : no_update_nodes) {
		index_U[alpha].insert(it_end, x);
		if (!beta_pointer) {
			for (int i = start_beta + 1; i <= beta; i++) pointer_U[alpha][i] = prev(it_end);
			beta_pointer = true;
		}
	}
	for (auto x : update_nodes) {
		index_U[alpha].insert(it_end, x);
		if (!beta_pointer) {
			for (int i = start_beta + 1; i <= beta; i++) pointer_U[alpha][i] = prev(it_end);
			beta_pointer = true;
		}
		if (!beta_plus_pointer) pointer_U[alpha][beta + 1] = prev(it_end), beta_plus_pointer = true;
	}
}
void BD_Index::inc_V(int beta, int alpha, unordered_set<int>& update_nodes) {
	if (update_nodes.empty()) return;
	if (alpha == beta) {
		for (auto x : update_nodes)
			index_V[beta].push_front(x);
		pointer_V[beta][alpha + 1] = index_V[beta].begin();
		return;
	}
	list<int>::iterator it_start = pointer_V[beta][alpha];
	check(pointer_V[beta].size() > alpha + 1, "pointer_V[beta].size() error");
	list<int>::iterator it_end = pointer_V[beta][alpha + 1];
	int start_alpha = alpha; while (start_alpha >= 0 && pointer_V[beta][start_alpha] == pointer_V[beta][alpha]) start_alpha--;

	unordered_set<int> no_update_nodes;
	while (it_start != it_end) {
		if (!update_nodes.count(*it_start))
			no_update_nodes.insert(*it_start);
		it_start = index_V[beta].erase(it_start);
	}
	bool alpha_pointer = false, alpha_plus_pointer = false;
	for (auto x : no_update_nodes) {
		index_V[beta].insert(it_end, x);
		if (!alpha_pointer) {
			for (int i = start_alpha + 1; i <= alpha; i++) pointer_V[beta][i] = prev(it_end);
			alpha_pointer = true;
		}
	}
	for (auto x : update_nodes) {
		index_V[beta].insert(it_end, x);
		if (!alpha_pointer) {
			for (int i = start_alpha + 1; i <= alpha; i++) pointer_V[beta][i] = prev(it_end);
			alpha_pointer = true;
		}
		if (!alpha_plus_pointer) pointer_V[beta][alpha + 1] = prev(it_end), alpha_plus_pointer = true;
	}
}
void BD_Index::dec_U(int alpha, int beta, unordered_set<int>& update_nodes) {
	if (update_nodes.empty()) return;
	if (beta == alpha) {
		int cnt = 0, len = update_nodes.size();
		for (auto it = pointer_U[alpha][beta]; cnt != len;) {
			if (update_nodes.count(*it))
				it = index_U[alpha].erase(it), cnt++;
			else
				it++;
		}
		pointer_U[alpha][beta] = index_U[alpha].begin();
		return;
	}
	list<int>::iterator it = pointer_U[alpha][beta];
	int start_beta = beta; while (start_beta >= 0 && pointer_U[alpha][start_beta] == pointer_U[alpha][beta]) start_beta--;
	list<int>::iterator it_end;
	if (beta != pointer_U[alpha].size() - 1) it_end = pointer_U[alpha][beta + 1];
	else it_end = index_U[alpha].end();
	while (it != it_end) {
		if (!update_nodes.count(*it))
			it++;
		else {
			if (it == pointer_U[alpha][beta]) {
				pointer_U[alpha][beta]++;
				it++;
			}
			else {
				index_U[alpha].insert(pointer_U[alpha][beta], *it);
				it = index_U[alpha].erase(it);
			}
		}
	}
	if (start_beta != beta - 1) {
		it = pointer_U[alpha][beta];
		int len = update_nodes.size();
		for (int i = 0; i < len; i++) it--;
		for (int i = start_beta + 1; i < beta; i++) pointer_U[alpha][i] = it;
	}
}
void BD_Index::dec_V(int beta, int alpha, unordered_set<int>& update_nodes) {
	if (update_nodes.empty()) return;
	if (alpha == beta + 1) {
		int cnt = 0, len = update_nodes.size();
		for (auto it = pointer_V[beta][alpha]; cnt != len;) {
			if (update_nodes.count(*it))
				it = index_V[beta].erase(it), cnt++;
			else
				it++;
		}
		pointer_V[beta][alpha] = index_V[beta].begin();
		return;
	}
	list<int>::iterator it = pointer_V[beta][alpha];
	int start_alpha = alpha; while (start_alpha >= 0 && pointer_V[beta][start_alpha] == pointer_V[beta][alpha]) start_alpha--;
	list<int>::iterator it_end;
	if (alpha != pointer_V[beta].size() - 1) it_end = pointer_V[beta][alpha + 1];
	else it_end = index_V[beta].end();
	while (it != it_end) {
		if (!update_nodes.count(*it))
			it++;
		else {
			if (it == pointer_V[beta][alpha]) {
				pointer_V[beta][alpha]++;
				it++;
			}
			else {
				index_V[beta].insert(pointer_V[beta][alpha], *it);
				it = index_V[beta].erase(it);
			}
		}
	}
	if (start_alpha != alpha - 1) {
		it = pointer_V[beta][alpha];
		int len = update_nodes.size();
		for (int i = 0; i < len; i++) it--;
		for (int i = start_alpha + 1; i < alpha; i++) pointer_V[beta][i] = it;
	}
}
void BD_Index::set_U(int alpha, int u) {
	unordered_map<int, int> neighbor;
	for (int j = 0; j < G.undeg[u]; j++) {
		Edge& ne = G.e[G.adj[u][j]];
		if (ne.deleted) continue;
		neighbor[ne.v] = alpha - 1;
	}
	int len = pointer_U[alpha].size();
	int now_r = alpha - 1;
	list<int>::iterator u_iterator;
	for (int b = alpha; b < len; b++) {
		list<int>::iterator it_start = pointer_U[alpha][b];
		list<int>::iterator it_end = b == len - 1 ? index_U[alpha].end() : pointer_U[alpha][b + 1];
		while (it_start != it_end) {
			if (neighbor.count(*it_start)) neighbor[*it_start] = b;
			if (*it_start == u) now_r = b, u_iterator = it_start;
			it_start++;
		}
	}
	vector<int> sorted_neighbor;
	for (auto p : neighbor) sorted_neighbor.push_back(p.second);
	sort(sorted_neighbor.begin(), sorted_neighbor.end(), greater<int>());
	int target_r = sorted_neighbor.size() > alpha ? sorted_neighbor[alpha] : alpha - 1;
	if (now_r == target_r) return;

	if (now_r > alpha - 1) {
		int start_beta = now_r;
		if (*(pointer_U[alpha][now_r]) == u) {
			while (start_beta >= alpha && *(pointer_U[alpha][start_beta]) == u) start_beta--;
			list<int>::iterator next_iterator = u_iterator; next_iterator++;
			for (int i = now_r; i > start_beta; i--) pointer_U[alpha][i] = next_iterator;
		}
		index_U[alpha].erase(u_iterator);
	}

	if (target_r > alpha - 1) {
		list<int>::iterator plus_one_iterator = target_r == pointer_U[alpha].size() - 1 ? index_U[alpha].end() : pointer_U[alpha][target_r + 1];
		index_U[alpha].insert(plus_one_iterator, u);
	}
}
void BD_Index::set_V(int beta, int v) {
	unordered_map<int, int> neighbor;
	for (int j = 0; j < G.undeg[v]; j++) {
		Edge& ne = G.e[G.adj[v][j]];
		if (ne.deleted) continue;
		neighbor[ne.u] = beta;
	}
	int len = pointer_V[beta].size();
	int now_r = beta - 1;
	list<int>::iterator v_iterator;
	for (int a = beta + 1; a < len; a++) {
		list<int>::iterator it_start = pointer_V[beta][a];
		list<int>::iterator it_end = a == len - 1 ? index_V[beta].end() : pointer_V[beta][a + 1];
		while (it_start != it_end) {
			if (neighbor.count(*it_start)) neighbor[*it_start] = a;
			if (*it_start == v) now_r = a, v_iterator = it_start;
			it_start++;
		}
	}
	vector<int> sorted_neighbor;
	for (auto p : neighbor) sorted_neighbor.push_back(p.second);
	sort(sorted_neighbor.begin(), sorted_neighbor.end(), greater<int>());
	int target_r = sorted_neighbor.size() > beta ? sorted_neighbor[beta] : beta;
	if (now_r == target_r) return;

	if (now_r > beta) {
		int start_alpha = now_r;
		if (*(pointer_V[beta][now_r]) == v) {
			while (start_alpha >= beta + 1 && *(pointer_V[beta][start_alpha]) == v) start_alpha--;
			list<int>::iterator next_iterator = v_iterator; next_iterator++;
			for (int i = now_r; i > start_alpha; i--) pointer_V[beta][i] = next_iterator;
		}
		index_V[beta].erase(v_iterator);
	}

	if (target_r > beta) {
		list<int>::iterator plus_one_iterator = target_r == pointer_V[beta].size() - 1 ? index_V[beta].end() : pointer_V[beta][target_r + 1];
		index_V[beta].insert(plus_one_iterator, v);
	}
}
void BD_Index::set_U_time(int alpha, int u) {
	int now_r = r_U[alpha][u];
	vector<int> sorted_neighbor;
	for (int j = 0; j < G.undeg[u]; j++) {
		if (G.e[G.adj[u][j]].deleted) continue;
		sorted_neighbor.push_back(r_U[alpha][G.e[G.adj[u][j]].v]);
	}
	sort(sorted_neighbor.begin(), sorted_neighbor.end(), greater<int>());
	int target_r = sorted_neighbor.size() > alpha ? sorted_neighbor[alpha] : -1;
	if (now_r == target_r) return;

	r_U[alpha][u] = target_r;
	if (now_r > alpha - 1) {
		list<int>::iterator u_iterator;
		list<int>::iterator it_start = pointer_U[alpha][now_r];
		list<int>::iterator it_end = now_r == pointer_U[alpha].size() - 1 ? index_U[alpha].end() : pointer_U[alpha][now_r + 1];
		while (it_start != it_end) {
			if (*it_start == u) {
				u_iterator = it_start;
				break;
			}
			it_start++;
		}
		int start_beta = now_r;
		if (*(pointer_U[alpha][now_r]) == u) {
			while (start_beta >= alpha && *(pointer_U[alpha][start_beta]) == u) start_beta--;
			list<int>::iterator next_iterator = u_iterator; next_iterator++;
			for (int i = now_r; i > start_beta; i--) pointer_U[alpha][i] = next_iterator;
		}
		index_U[alpha].erase(u_iterator);
	}

	if (target_r > alpha - 1) {
		list<int>::iterator plus_one_iterator = target_r == pointer_U[alpha].size() - 1 ? index_U[alpha].end() : pointer_U[alpha][target_r + 1];
		index_U[alpha].insert(plus_one_iterator, u);
	}
}
void BD_Index::set_V_time(int beta, int v) {
	int now_r = r_V[beta][v];
	vector<int> sorted_neighbor;
	for (int j = 0; j < G.undeg[v]; j++) {
		if (G.e[G.adj[v][j]].deleted) continue;
		sorted_neighbor.push_back(r_V[beta][G.e[G.adj[v][j]].u]);
	}
	sort(sorted_neighbor.begin(), sorted_neighbor.end(), greater<int>());
	int target_r = sorted_neighbor.size() > beta ? sorted_neighbor[beta] : -1;
	if (now_r == target_r) return;

	r_V[beta][v] = target_r;
	if (now_r > beta) {
		list<int>::iterator v_iterator;
		list<int>::iterator it_start = pointer_V[beta][now_r];
		list<int>::iterator it_end = now_r == pointer_V[beta].size() - 1 ? index_V[beta].end() : pointer_V[beta][now_r + 1];
		while (it_start != it_end) {
			if (*it_start == v) {
				v_iterator = it_start;
				break;
			}
			it_start++;
		}
		int start_alpha = now_r;
		if (*(pointer_V[beta][now_r]) == v) {
			while (start_alpha >= beta + 1 && *(pointer_V[beta][start_alpha]) == v) start_alpha--;
			list<int>::iterator next_iterator = v_iterator; next_iterator++;
			for (int i = now_r; i > start_alpha; i--) pointer_V[beta][i] = next_iterator;
		}
		index_V[beta].erase(v_iterator);
	}

	if (target_r > beta) {
		list<int>::iterator plus_one_iterator = target_r == pointer_V[beta].size() - 1 ? index_V[beta].end() : pointer_V[beta][target_r + 1];
		index_V[beta].insert(plus_one_iterator, v);
	}
}
void BD_Index::compute_memory() {
	size_t sum = 0;

	sum += sizeof(vector<vector<list<int>::iterator> >);
	sum += pointer_U.capacity() * sizeof(vector<list<int>::iterator>);
	for (auto x : pointer_U)
		sum += x.capacity() * sizeof(list<int>::iterator);
	sum += sizeof(vector<list<int> >);
	sum += index_U.capacity() * sizeof(list<int>);
	for (auto x : index_U)
		sum += x.size() * (sizeof(int) + 2 * sizeof(void*));

	sum += sizeof(vector<vector<list<int>::iterator> >);
	sum += pointer_U.capacity() * sizeof(vector<list<int>::iterator>);
	for (auto x : pointer_U)
		sum += x.capacity() * sizeof(list<int>::iterator);
	sum += sizeof(vector<list<int> >);
	sum += index_U.capacity() * sizeof(list<int>);
	for (auto x : index_U)
		sum += x.size() * (sizeof(int) + 2 * sizeof(void*));

	printf("- %-20s: %.4lf\n", "BD_Index memory (MB)", 1.0 * sum / 1024 / 1024);
	// printf("- %-20s: %.4lf\n", "BD_Index memory (GB)", 1.0 * sum / 1024 / 1024 / 1024);

	if (algorithm_used == TIME) {
		sum = 0;
		sum += size_t(6) * (G.pseudoarboricity + 1) * G.N * sizeof(int);
		sum += size_t(2) * (G.pseudoarboricity + 1) * G.M * sizeof(bool);
		printf("- %-20s: %.4lf\n", "Auxiliary memory (MB)", 1.0 * sum / 1024 / 1024);
		// printf("- %-20s: %.4lf\n", "Auxiliary memory (GB)", 1.0 * sum / 1024 / 1024 / 1024);
	}
}
void BD_Index::display() {
	vector<int> now_log;
	bool output = false;
	if (output) printf("--------------------\n");
	for (int alpha = 0; alpha <= G.pseudoarboricity; alpha++) {
		for (int beta = alpha; beta < pointer_U[alpha].size(); beta++) {
			int num;
			if (beta != pointer_U[alpha].size() - 1)
				num = distance(pointer_U[alpha][beta], pointer_U[alpha][beta + 1]);
			else
				num = distance(pointer_U[alpha][beta], index_U[alpha].end());
			if (output) printf("D(%d, %d) = %d, ", alpha, beta, num);
			now_log.push_back(num);
		}
		if (output) printf("\n");
	}
	for (int beta = 0; beta <= G.pseudoarboricity; beta++) {
		for (int alpha = beta + 1; alpha < pointer_V[beta].size(); alpha++) {
			int num;
			if (alpha != pointer_V[beta].size() - 1)
				num = distance(pointer_V[beta][alpha], pointer_V[beta][alpha + 1]);
			else
				num = distance(pointer_V[beta][alpha], index_V[beta].end());
			if (output) printf("D(%d, %d) = %d, ", alpha, beta, num);
			now_log.push_back(num);
		}
		if (output) printf("\n");
	}
	log.push_back(now_log);
}
void BD_Index::check_correctness() {
	int len = log.size();
	for (int i = 0; i < len / 2; i++) {
		check(log[i].size() == log[len - 1 - i].size(), "log size error");
		for (int j = 0; j < log[i].size(); j++)
			check(log[i][j] == log[len - 1 - i][j], "log error");
	}
	return;
}
void BD_Index::solve_query(int alpha, int beta) {
	answer.clear();
	if (beta >= alpha) {
		if (alpha < pointer_U.size()) {
			if (beta < pointer_U[alpha].size()) {
				auto start = pointer_U[alpha][beta];
				while (start != index_U[alpha].end()) {
					answer.push_back(*start);
					start++;
				}
			}
		}
	}
	else {
		if (beta < pointer_V.size()) {
			if (alpha < pointer_V[beta].size()) {
				auto start = pointer_V[beta][alpha];
				while (start != index_V[beta].end()) {
					answer.push_back(*start);
					start++;
				}
			}
		}
	}
}
void analyze_instruction_timings() {
	int n = instructions.size();
	check(n == (int)process_time.size(), "Error: instructions.size() != process_time.size()\n");

	vector<double> waiting_times;
	waiting_times.reserve(n);

	double current_time = 0.0;

	for (int i = 0; i < n; i++) {
		const auto& inst = instructions[i];

		double arrival = double(inst.time) / 1000;
		double service = process_time[i];

		if (current_time < arrival) {
			current_time = arrival;
		}

		current_time += service;

		double wait = current_time - arrival;
		waiting_times.push_back(wait);
	}

	double query_sum = 0.0, query_max = 0.0, update_sum = 0.0, update_max = 0.0;
	int query_num = 0, update_num = 0;
	for (int i = 0; i < n; i++) {
		if (instructions[i].instruction_type == QUERY) {
			query_num++;
			query_sum += waiting_times[i];
			query_max = max(query_max, waiting_times[i]);
		}
		else if (instructions[i].instruction_type == INSERT || instructions[i].instruction_type == DELETE) {
			update_num++;
			update_sum += waiting_times[i];
			update_max = max(update_max, waiting_times[i]);
		}
	}
	double query_avg = query_num == 0 ? 0.0 : query_sum / query_num;
	double update_avg = update_num == 0 ? 0.0 : update_sum / update_num;

	cout << fixed << setprecision(6);
	cout << "Analyzed " << waiting_times.size() << " UPDATE/QUERY instructions\n";
	cout << "Average query turnaround time: " << query_avg << " s\n";
	cout << "Max query turnaround time:     " << query_max << " s\n";
	cout << "Average update turnaround time: " << update_avg << " s\n";
	cout << "Max update turnaround time:     " << update_max << " s\n";
}


int main(int argc, char** argv) {
	if (argc != 4) {
	argument_error:
		printf("Usage: ./main <dataset_address> <update_algorithm> <instruction_address>\n");
		printf("update algorithm:\n");
		printf("-space: space efficient\n");
		printf("-time: time efficient\n");
		return 0;
	}
	Timer timer; double runtime;
	char dataset_address[1000]; strcpy(dataset_address, argv[1]);
	if (strcmp(argv[2], "-space") == 0) algorithm_used = SPACE;
	else if (strcmp(argv[2], "-time") == 0) algorithm_used = TIME;
	else goto argument_error;

	printf("----------Now processing %s----------\n", dataset_address);
	printf("- %-20s: %s\n", "Algorithm used", argv[2]);

	timer.start();
	G.read_graph_from_dataset(dataset_address);
	timer.end(); runtime = timer.time();
	printf("- %-20s: %d, %d, %d\n", "|E|, |U|, |V|", G.M, G.U, G.V);
	printf("- %-20s: %lf\n", "Read graph time", runtime);

	timer.start();
	G.construct_index();
	timer.end(); runtime = timer.time();
	printf("- %-20s: %lf\n", "Construct bd time", runtime);
	printf("- %-20s: %d\n", "Pseudoarboricity", G.pseudoarboricity);

	bd.compute_memory();

	for (int alpha = 0; alpha < bd.index_U.size(); alpha++) {
		bd.index_U[alpha].clear();
		for (int beta = alpha; beta < bd.pointer_U[alpha].size(); beta++) {
			bd.pointer_U[alpha][beta] = bd.index_U[alpha].end();
		}
	}
	for (int beta = 0; beta < bd.index_V.size(); beta++) {
		bd.index_V[beta].clear();
		for (int alpha = beta + 1; alpha < bd.pointer_V[beta].size(); alpha++) {
			bd.pointer_V[beta][alpha] = bd.index_V[beta].end();
		}
	}
	bd.log.clear();
	int len = bd.r_U.size();
	for (int p = 0; p < len; p++) {
		for (int x = 0; x < G.N; x++) {
			bd.indeg_U[p][x] = bd.neighbor_cnt_U[p][x] = bd.indeg_V[p][x] = bd.neighbor_cnt_V[p][x] = 0;
			bd.r_U[p][x] = bd.r_V[p][x] = -1;
		}
	}
	for (int i = 0; i < G.M; i++) G.e[i].deleted = true;

	char instruction_address[1000]; strcpy(instruction_address, argv[3]);
	get_instructions(instruction_address);

	int instruction_number = instructions.size();
	for (int i = 0; i < instruction_number; i++) {
        // if (algorithm_used == SPACE && i % 1000 == 0) printf("i = %d\n", i);
		timer.start();
		Instruction& now_instruction = instructions[i];
		if (now_instruction.instruction_type == INSERT) {
			if (algorithm_used == SPACE) G.insert_edge_space(now_instruction.t1);
			else if (algorithm_used == TIME) G.insert_edge_time(now_instruction.t1);
		}
		else if (now_instruction.instruction_type == DELETE) {
			if (algorithm_used == SPACE) G.delete_edge_space(now_instruction.t1);
			else if (algorithm_used == TIME) G.delete_edge_time(now_instruction.t1);
		}
		else if (now_instruction.instruction_type == QUERY) {
			bd.solve_query(now_instruction.t1, now_instruction.t2);
		}
		timer.end();
		process_time[i] = timer.time();
	}

	analyze_instruction_timings();

	return 0;
}