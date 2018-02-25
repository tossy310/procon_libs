// 2部グラフのマッチング O(VE), Dinicより遅め
class BipartiteMatching {
private:
  int n;
  vector<vector<int>> graph;
  vector<int> match;
  vector<bool> used;
  bool dfs(int v){
    used[v] = true;
    for(int u : graph[v]){
      int w = match[u];
      if(w<0 || (!used[w] && dfs(w))){
        match[v] = u;
        match[u] = v;
        return true;
      }
    }
    return false;
  }
public:
  BipartiteMatching(int _n) : n(_n){
    graph.resize(n);
    match.resize(n);
    used.resize(n);
  }
  void add_edge(int u, int v){ // 0-indexed
    graph[u].pb(v);
    graph[v].pb(u);
  }
  int matching(){ // matching size
    int res = 0;
    fill(begin(match), end(match), -1);
    rep(v, n){
      if(match[v] < 0){
        fill(begin(used), end(used), false);
        if(dfs(v)) res++;
      }
    }
    return res;
  }
};

/*
  2部グラフの場合，
    |最大マッチング| = |最小点カバー|
  一般に，
    |最大安定集合| + |最小点カバー| = |V|
  また，各連結集合に対して，
    |最大マッチング| + |最小辺カバー| = |V|
*/


// 2部マッチングの片方が[0,left_n)のときの最小頂点被覆
// Dinicのときは (sから到達不可な左側頂点) ⋃ (sから到達可能な右側頂点) が答えなので，残余グラフをsからdfsで求まる．
  vector<bool> color;
  void color_dfs(int v, int left_n){
    color[v] = true;
    for(auto to : graph[v]){
      if((v < left_n)^(match[v]==to)){
        if(!color[to]) color_dfs(to, left_n);
      }
    }
  }
  vector<int> minEdgeCover(int left_n){
    color.resize(n,false);
    rep(i,left_n) if(match[i]==-1) color_dfs(i, left_n);
    vector<int> res;
    rep(i,n) if((i<left_n) ^ color[i]) res.pb(i);
    return res;
  }
