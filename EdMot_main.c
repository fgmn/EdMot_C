
#include <stdio.h>
#include <string.h>
#include <math.h>

#define N 120
#define max(a, b) ((a) >= (b) ? (a) : (b))
#define min(a, b) ((a) <= (b) ? (a) : (b))

/************ �������� ***************/
int A[N][N], W[N][N];	// �ڽӾ���(�Գ�)��ͼ��WΪģ�����
int n, m;		// �ڵ���������
int vis[N];		// �ڵ��Ƿ���ʹ�
int com[N];		// �ڵ�������ͨ����
int com_sz[N];	// ��ͨ��Ĵ�С���������ڵ�����

// Louvain
struct node { int a, b; } Node[N];
//int a_new[N][N] = { 0 }, a_old[N][N] = { 0 };
int COM[N];	// todo
//int ind_com_now[N], ind_com_full[N];
//int cnt_now[N], cnt_full[N];
//int e;

int std_label[N];
/************ �������� ***************/

/************ �������� ***************/
// ����
void input();

// �������ͨ����
void find_LCC(int(*a)[N]);
void dfs(int u, int c, int(*a)[N]);
// ��ģ�����
void get_Motif_mat();

// partition
void partition(int(*a)[N], int n_a, int k, char* meth);
void kmeans_Cluster(int(*a)[N], int k);
void Louvain(int(*a)[N], int n_a, int ctrl);
double compute_modularity(int(*a)[N], int n_a, int* label);
void reindex_com(int* com_old);

// �����ǿ��
void strengthen_edge();

// ָ��
void NMI(int *std_label, int *cluster_label);
double Modulatity(int(*a)[N], int* label);

/************ �������� ***************/

int main() {
	printf("EdMot project:\n");

	input();

	// ��Ԥ������ԭ������������ͨ������ֱ���ų�ԭͼ�еĹ����㣬���ǿ϶���������������γ�
	find_LCC(A);
	//for (int i = 1; i <= n; i++) {
	//	for (int j = 1; j <= n; j++) {
	//		printf("%d ", A[i][j]);
	//	}
	//	printf("\n");
	//}
	printf("\nLouvain\n");
	partition(A, n, 0, "Louvain");
	NMI(std_label, COM);

	get_Motif_mat();
	//for (int i = 1; i <= n; i++) {
	//	for (int j = 1; j <= n; j++) {
	//		printf("%3d", W[i][j]);
	//	}
	//	printf("\n");
	//}

	// ��ͼ��������ͨ����
	find_LCC(W);
	//for (int i = 1; i <= n; i++) {
	//	for (int j = 1; j <= n; j++) {
	//		if (W[i][j] > 1) W[i][j] = 1;
	//	}
	//}

	// ��һ�� partition
	printf("\nMotif-Louvain\n");
	partition(W, n, 0, "Louvain");
	// ����ֵ�洢�� COM
	NMI(std_label, COM);
	// �����ǿ�߼��ϣ�����Եǿ��
	strengthen_edge();

	// �ڶ��� partition
	printf("\nEdMot-Louvain\n");
	partition(A, n, 0, "Louvain");
	NMI(std_label, COM);

	return 0;
}

void input() {
	FILE* fp = NULL;
	int err = fopen_s(&fp, "./polbooks.txt", "r");
	if (err) { printf("%d", err); return; }
	fscanf_s(fp, "%d %d", &n, &m);
	//printf("%d %d\n", n, m);

	int x, y;
	for (int i = 1; i <= m; i++) {
		fscanf_s(fp, "%d %d", &x, &y);
		//printf("%d %d\n", x, y);
		A[x][y] = A[y][x] = 1;
	}

	// ��ȡ��ǩ
	err = fopen_s(&fp, "./label.txt", "r");
	if (err) { printf("%d", err); return; }
	for (int i = 1; i <= n; i++) {
		fscanf_s(fp, "%d", &x);
		std_label[i] = x;
	}
}


void find_LCC(int(*a)[N]) {
	int cnt_com = 0;
	memset(vis, 0, sizeof vis);
	for (int i = 1; i <= n; i++) {
		if (!vis[i]) dfs(i, ++cnt_com, a);
		//printf("%d ", com[i]);
	}
	//printf("\n");

	// ���������ͨ������Ӧ�ı��
	int com_sz_max = 0, com_num = 0;
	for (int i = 1; i <= cnt_com; i++) {
		if (com_sz[i] > com_sz_max) {
			com_num = i;
			com_sz_max = com_sz[i];
		}
	}
	// �����������ͨ����
	for (int i = 1; i <= n; i++) {
		for (int j = i + 1; j <= n; j++) {
			if (com[i] != com_num && com[j] != com_num) {
				a[i][j] = a[j][i] = 0;
			}
		}
	}
}

void dfs(int u, int c, int(*a)[N]) {
	vis[u] = 1;
	com[u] = c;
	com_sz[c]++;
	for (int i = 1; i <= n; i++) {
		if (!vis[i] && A[u][i]) dfs(i, c, a);
	}
}

// ����ģ�����ѡ������ģ�壬����Ϊm4�������http://snap.stanford.edu/higher-order/code.html��
void get_Motif_mat() {
	for (int i = 1; i <= n; i++) {
		for (int j = i + 1; j <= n; j++) {
			for (int k = j + 1; k <= n; k++) {
				if (A[i][j] + A[j][k] + A[k][i] == 3) {
					W[i][j]++; W[j][i]++;
					W[j][k]++; W[k][j]++;
					W[k][i]++; W[i][k]++;
				}
			}
		}
	}
}

void partition(int(*a)[N], int n_a, int k, char* meth) {
	if (!strcmp(meth, "Louvain")) {
		//printf("Louvain\n");
		Louvain(a, n_a, 0);
	}
	if (!strcmp(meth, "SC")) {
		kmeans_Cluster(a, k);
	}
	// ���� partition ����
}

void kmeans_Cluster(int(*a)[N], int k) {}

// ���������㷨 Louvain���ɲο���https://blog.csdn.net/qq_16543881/article/details/122825957��
// matlab ����ֻ������ģ����Ż��׶�
void Louvain(int(*a)[N], int n_a, int ctrl) {
	int Niter = 1;	// �����ִ�
	double mat_sum = 0;
	double mat_rowsum[N] = { 0 };
	for (int i = 1; i <= n_a; i++) {
		for (int j = 1; j <= n_a; j++) {
			mat_rowsum[i] += a[i][j];
		}
		mat_sum += mat_rowsum[i];
	}
	double K[N] = { 0 }, sum_tot[N] = { 0 }, sum_in[N] = { 0 };
	// K���������ڵ��Ȩ���ܺ� sum_tot������������Ȩ�غ� sum_in�������ڲ���Ȩ�غ�
	int com[N] = { 0 };
	// ��ʼ����ÿ���ڵ���Ϊһ������
	for (int i = 1; i <= n_a; i++) {
		sum_tot[i] = K[i] = mat_rowsum[i];
		sum_in[i] = a[i][i];
		com[i] = i;
	}

	double sum_cost = 10;
	int gain = 1;

	while (gain) {
		double cost[N] = { 0 };
		gain = 0;
		for (int i = 1; i <= n_a; i++) {
			int ci = com[i], c_new = ci, c_new_tmp = 0;
			com[i] = -1;
			double gain_vec[N] = { 0 };	// ���ڵ� i �ϲ�������������
			double best_inc = -1;		// ���ŵĺϲ�����

			// ���ڵ� i �Ƴ�ԭ���� ci
			sum_tot[ci] -= K[i];
			// ��ȡͬһ�����������ڵ�
			int cni[N] = { 0 }, cnt_cni = 0;
			double sum_cn = 0;
			for (int j = 1; j <= n_a; j++) {
				if (com[j] == ci) { 
					cni[++cnt_cni] = j;
					sum_cn += a[i][j];
				}
			}
			// �˴���֤ ci <= n_a
			sum_in[ci] -= 2 * sum_cn + a[i][i];

			// �����ڽӵ㣬���Ժϲ�
			for (int j = 1; j <= n_a; j++) {
				if (a[i][j]) {
					int cj = com[j];
					if (gain_vec[j] < 1e-10) {
						int cnj[N] = { 0 }, cnt_cnj = 0;
						double ki_in = 0;	// �ڵ� i ������ cj �Ĺ���Ȩ��
						for (int k = 1; k <= n_a; k++) {
							if (com[k] == cj) {
								cnj[++cnt_cnj] = k;
								ki_in += a[i][k];
							}
						}
						// ����
						gain_vec[cj] = 2 * ki_in / mat_sum - 2 * K[i] * sum_tot[cj] / (mat_sum * mat_sum);
					}
					if (gain_vec[cj] > best_inc) {
						best_inc = gain_vec[cj];
						c_new_tmp = cj;	// ���ŵĺϲ�ѡ��
					}
				}
			}
			if (best_inc > 0) {
				cost[i] = best_inc;
				c_new = c_new_tmp;
			}
			// ������ѡ����кϲ�
			int ck[N] = { 0 }, cnt_ck = 0;
			double sum_ck = 0;
			for (int j = 1; j <= n_a; j++) {
				if (com[j] == c_new) {
					ck[++cnt_ck] = j;
					sum_ck += a[i][j];
				}
			}
			sum_in[c_new] += 2 * sum_ck;
			sum_tot[c_new] += K[i];
			com[i] = c_new;

			if (c_new != ci) { gain = 1; }
		}

		for (int i = 1; i <= n_a; i++) {
			sum_cost += cost[i];
		}
		reindex_com(com);
		double mod = compute_modularity(a, n_a, com);
		//printf("It %d - Mod=%f\n", Niter, mod);
		Niter++;
	}

	//printf("\n");
	printf("label as follows:\n");
	for (int i = 1; i <= n_a; i++) {
		COM[i] = com[i];
		printf("%d ", com[i]);
	}
	printf("\n");

	/*
	Niter--;
	
	if (ctrl) {
		for (int i = 1; i <= n_a; i++) {
			for (int j = 1; j <= n_a; j++) {
				a_new[i][j] = a_old[i][j] = a[i][j];
			}
		}
		int com_now[N] = { 0 }, com_full[N] = { 0 };
		for (int i = 1; i <= n; i++) com_now[i] = com_full[i] = com[i];
		int k = 2;

		while (1) {
			int n_node = 0, n = 0;	// todo
			for (int i = 1; i <= n_a; i++) {
				for (int j = 1; j <= n_a; j++) {
					a_old[i][j] = a_new[i][j];
				}
			}
			int cnt[N] = { 0 }, n_com = 0;	// ͳ��������
			for (int i = 1; i <= n; i++) {
				if (!cnt[com_now[i]]) n_com++;
				cnt[com_now[i]]++;
			}
			for (int i = 1; i <= n_com; i++) {
				for (int j = 1; j <= n_node; j++) {
					if (com_now[j] == i) {
						ind_com_now[++cnt_now[i]] = j;
					}
				}
			}
			for (int i = 1; i <= n_com; i++) {
				for (int j = 1; j <= n; j++) {
					if (com_full[j] == i) {
						ind_com_full[++cnt_full[i]] = j;
					}
				}
			}
			// ��������Ϊһ���ڵ�
			for (int i = 1; i <= n_com; i++) {
				for (int j = i; j <= n_com; j++) {
					double sum_ij = 0;
					for (int p = 1; p <= cnt_now[i]; p++) {
						for (int q = 1; q <= cnt_now[j]; q++) {
							sum_ij += a_old[ind_com_now[p]][ind_com_now[q]];
						}
					}
					a_new[i][j] = a_new[j][i] = sum_ij;
				}
			}
			Louvain(a_new, n_com, 0);    
			if (e != 1) {
				for (int i = 1; i <= n_a; i++) {
					com_full[i] = 0;
					com_now[i] = com[i];	// com �洢����ֵ
				}
				for (int i = 1; i <= n_com; i++) {
					for (int j = 1; j <= cnt_full[i]; j++) {
						com_full[ind_com_full[j]] = com_now[i];
					}
				}
				reindex_com(com_full);
				int is_same = 1;
				for (int i = 1; i <= n_a; i++) {
					if (COM[k - 1][i] != com_full[i]) {
						is_same = 0;
						break;
					}
				}
				if (is_same) { return; }
			}
			else {
				return;
			}
			k++;
		}
	}*/
}

// ����ģ��ȣ��ж��ּ��㷽ʽ��https://zhuanlan.zhihu.com/p/178790546��
// ���Ĺ�ʽ�� matlab ʵ�ֲ�ͬ����ѭ����
double compute_modularity(int(*a)[N], int n_a, int *label) {
	double mod = 0, mat_sum = 0;
	double mat_rowsum[N] = { 0 };	// �������ÿһ�еĺͣ����ڼ���
	for (int i = 1; i <= n_a; i++) {
		for (int j = 1; j <= n_a; j++) {
			mat_rowsum[i] += a[i][j];
		}
		mat_sum += mat_rowsum[i];
	}
	
	int com_num = 0;	// ������
	for (int i = 1; i <= n_a; i++) {
		com_num = max(com_num, label[i]);
	}
	for (int i = 1; i <= com_num; i++) {
		int Ci[N] = { 0 }, cnt = 0;
		// ��ȡ�������еĽڵ�
		for (int j = 1; j <= n; j++) {
			if (label[j] == i) {
				Ci[++cnt] = j;
			}
		}
		double Ec = 0, Et = 0;
		for (int j = 1; j <= cnt; j++) {
			Et += mat_rowsum[Ci[j]];
		}
		for (int j = 1; j <= cnt; j++) {
			for (int k = 1; k <= cnt; k++) {
				Ec += a[Ci[j]][Ci[k]];
			}
		}
		if (Et > 0) {
			mod += Ec / mat_sum - (Et / mat_sum) * (Et / mat_sum);
		}
	}
	return mod;
} 

// �ϲ�֮�󣬲������ȱʧ��������±���
void reindex_com(int *com_old) {
	int cnt[N] = { 0 }, cnt_unique = 0;
	for (int i = 1; i <= n; i++) {
		cnt[com_old[i]]++;
	}
	for (int i = 1; i <= n; i++) {
		if (cnt[i]) cnt_unique++;
		Node[i].a = cnt[i];
		Node[i].b = i;
	}
	int com_new[N] = { 0 };
	// ð������
	for (int i = 1; i <= n; i++) {
		for (int j = i + 1; j <= n; j++) {
			if (Node[i].a < Node[j].a) {
				struct node tmp = Node[j];
				Node[j] = Node[i];
				Node[i] = tmp;
			}
		}
	}
	for (int i = 1; i <= cnt_unique; i++) {
		for (int j = 1; j <= n; j++) {
			if (com_old[j] == Node[i].b) {
				com_new[j] = i;
			}
		}
	}
	for (int i = 1; i <= n; i++) {
		com_old[i] = com_new[i];
	}
}

void strengthen_edge() {
	//printf("strength edges:\n");
	for (int i = 1; i <= n; i++) {
		for (int j = i + 1; j <= n; j++) {
			if (COM[i] == COM[j]) {
				//if (!A[i][j]) printf("(%d, %d) ", i, j);
				A[i][j] = A[j][i] = 1;
			}
		}
	}
	printf("\n");
}

// ����Ч������ָ�꣬��׼������Ϣ NMI ���ο���https://blog.csdn.net/qq_42122496/article/details/106193859��
void NMI(int* std_label, int* cluster_label) {
	int u1[N] = { 0 }, u2[N] = { 0 };
	int cnt_u1 = 0, cnt_u2 = 0;
	int cnt1[N] = { 0 }, cnt2[N] = { 0 };
	for (int i = 1; i <= n; i++) {
		cnt1[std_label[i]]++;
		cnt2[cluster_label[i]]++;
	}
	for (int i = 1; i <= n; i++) {
		if (cnt1[i] > 0) u1[++cnt_u1] = i;
		if (cnt2[i] > 0) u2[++cnt_u2] = i;
	}
	double I = 0;
	for (int i = 1; i <= cnt_u1; i++) {
		for (int j = 1; j <= cnt_u2; j++) {
			double nc = 0;
			for (int k = 1; k <= n; k++) {
				if (std_label[k] == i && cluster_label[k] == j) { nc += 1; }
			}
			if (nc < 1e-2) continue;
			I += (nc / n) * (log(n * nc / cnt1[u1[i]] / cnt2[u2[j]]));
		}
	}
	double entropy1 = 0, entropy2 = 0;
	for (int i = 1; i <= cnt_u1; i++) {
		entropy1 -= (1.0 * cnt1[u1[i]] / n) * log(1.0 * cnt1[u1[i]] / n);
	}
	for (int i = 1; i <= cnt_u2; i++) {
		entropy2 -= (1.0 * cnt2[u2[i]] / n) * log(1.0 * cnt2[u2[i]] / n);
		//printf("1.0 * cnt2[u2[i]] / n = %lf\n", 1.0 * cnt2[u2[i]] / n);
	}

	double z = 2 * I / (entropy1 + entropy2);
	printf("NMI: %lf\n", z);
}

double Modulatity(int(*a)[N], int *label) {}