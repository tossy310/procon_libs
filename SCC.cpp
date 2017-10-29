// Strongly Connected Components O(E + V)
class SCC {
private:
  vector<vector<int>> &G;
  vector<vector<int>> &rG;
  vector<bool> used;
  vector<int> vs;
  int n;
  void dfs(int v){
    used[v] = true;
    for(auto to : G[v]) if(!used[to]) dfs(to);
    vs.pb(v);
  }
  void rdfs(int v, int k){
    used[v] = true;
    cmp[v] = k;
    group.back().pb(v);
    for(auto to : rG[v]) if(!used[to]) rdfs(to, k);
  }
public:
  vector<int> cmp; //頂点iが属するSCCのID
  vector<vector<int>> group; // SCCのトポロジカル順i番目に属する頂点番号群
  int gc = 0; // group count
  SCC(vector<vector<int>> &g, vector<vector<int>> &r) : G(g), rG(r){
    n = G.size();
    used = vector<bool>(n, false);
    cmp.resize(n);
    rep(i,n) if(!used[i]) dfs(i);
    fill(all(used), false);
    for(int i=vs.size()-1; i>=0; i--) if(!used[vs[i]]){
      group.pb(vector<int>());
      rdfs(vs[i], gc++);
    }
  }
  inline bool same(int i, int j){ return cmp[i]==cmp[j]; }
  inline int operator[](int i){ return cmp[i]; }
};

// 2-SAT O(M + N) M:# of vars, N:# of clauses
class TwoSat {
private:
  vector<vector<int> > vec;
  vector<vector<int> > rev;
  vector<bool> res;
  int v;
public:
  TwoSat(int n) : v(n){
    vec.resize(2*v);
    rev.resize(2*v);
    res.resize(v);
  }
  inline void add_edge(int a, bool ab, int b, bool bb){
    // add (a -> b)
    int va = a + (ab?0:v);
    int vb = b + (bb?0:v);
    vec[va].pb(vb);
    rev[vb].pb(va);
  }
  inline void add(int a, bool ab, int b, bool bb){
    // a or b <-> ((!a -> b) and (!b -> a))
    add_edge(a, !ab, b, bb);
    add_edge(b, !bb, a, ab);
  }
  bool exec(){
    SCC scc(vec, rev);
    rep(i,v){
      if(scc.same(i, i+v)) return false;
      res[i] = scc[i]>scc[i+v];
    }
    return true;
  }
  inline bool operator[](int i){ return res[i]; }
};


// for POJ
class SCC {
private:
public:
  vector<vector<int> > &G;
  vector<vector<int> > &rG;
  vector<bool> used;
  vector<int> vs;
  vector<int> cmp;
  int n;
  int gc; // group count
  void dfs(int v){
    used[v] = true;
    rep(i,G[v].size()) if(!used[G[v][i]]) dfs(G[v][i]);
    vs.pb(v);
  }
  void rdfs(int v, int k){
    used[v] = true;
    cmp[v] = k;
    rep(i,rG[v].size()) if(!used[rG[v][i]]) rdfs(rG[v][i], k);
  }
  SCC(vector<vector<int> > &g, vector<vector<int> > &r) : G(g), rG(r){
    n = G.size();
    used = vector<bool>(n,false);
    vs.clear();
    cmp.resize(n);
    rep(i,n) if(!used[i]) dfs(i);
    fill(all(used), false);
    gc=0;
    for(int i = vs.size()-1; i>=0; i--) if(!used[vs[i]]) rdfs(vs[i], gc++);
  }
  inline bool same(int i, int j){ return cmp[i]==cmp[j]; }
  inline int operator[](int i){ return cmp[i]; }
};
