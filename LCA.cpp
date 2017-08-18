vector<int> tree[100010];

class LCA {
public:
  int n,ln; // number of nodes and log n
  vector<vector<int> > parent;
  vector<int> depth;
  LCA(int _n, int root=-1) : n(_n), depth(_n){
    ln=0;
    while(n>(1<<ln)) ln++;  // calc log n
    parent = vector<vector<int> >(ln, vector<int>(n));
    if(root!=-1) init(root);
  }
  void dfs(int v, int p, int d){
    parent[0][v]=p;
    depth[v]=d;
    rep(i,tree[v].size()) if(tree[v][i]!=p) dfs(tree[v][i], v, d+1);
  }
  void init(int root){
    dfs(root, -1, 0);
    for(int k=0; k+1<ln; k++){
      for(int v=0; v<n; v++){
        if(parent[k][v] < 0) parent[k+1][v] = -1;
        else parent[k+1][v] = parent[k][parent[k][v]];
      }
    }
  }
  int query(int u, int v){
    if(depth[u] > depth[v]) swap(u,v);
    for(int k=0; k<ln; k++){
      if((depth[v]-depth[u])>>k & 1) v = parent[k][v];
    }
    if(u==v) return u;
    for(int k=ln-1; k>=0; k--){
      if(parent[k][u] != parent[k][v]){
        u = parent[k][u];
        v = parent[k][v];
      }
    }
    return parent[0][u];
  }
};
