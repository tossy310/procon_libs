// 2部グラフのマッチング
// O(VE)
class BipartiteMatching{
private:
  int n; //頂点数
  vector<vector<int> > graph; //グラフ
  vector<int> match; //マッチングペア
  vector<bool> used;
  bool dfs(int v){
    used[v]=true;
    rep(i, graph[v].size()){
      int u = graph[v][i], w = match[u];
      if(w<0 || (!used[w] && dfs(w)) ){ //増加路が存在したら，つなぎかえる
        match[v]=u;
        match[u]=v;
        return true;
      }
    }
    return false;
  }
public:
  BipartiteMatching(int nn) : n(nn){
    graph.resize(n);
    match.resize(n);
    used.resize(n);
  }
  void add_edge(int u, int v){ //0-indexedでマッチング辺の追加
    graph[u].pb(v);
    graph[v].pb(u);
  }
  int matching(){ //マッチング数を返す
    int res=0;
    fill(all(match), -1);
    rep(v, n){
      if(match[v]<0){
        fill(all(used), false);
        if(dfs(v)) res++; // 増加路が存在したので，マッチングが増加
      }
    }
    return res;
  }
}; //END class BipartiteMatching
