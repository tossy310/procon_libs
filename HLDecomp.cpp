// TODO 森の場合に複数の木を一気に管理．
// see http://beet-aizu.hatenablog.com/entry/2017/12/12/235950
// verified AOJ GRL 5C(LCA), 5D(for_each_edge), 5E(for_each_edge)

template<int SZ>
class HLDecomp {
private:
  int n;
  void dfs(const int root){
    stack<pair<int,bool> > st;
    par[root] = -1;
    dep[root] = 0;
    st.push({root,false});
    while(!st.empty()){
      int v = st.top().first;
      bool &b = st.top().second;
      if(!b){
        // initial visit of v
        b = true;
        for(int u : tree[v]) if(u != par[v]){
          par[u] = v;
          dep[u] = dep[v] + 1;
          st.push({u, false});
        }
      }
      else {
        // second visit
        st.pop();
        int cur_max = 0;
        for(int u : tree[v]) if(u != par[v]){
          sub_sz[v] += sub_sz[u];
          if(sub_sz[u] > cur_max){
            cur_max = sub_sz[u];
            heavy_child[v] = u;
          }
        }
      }
    }
  }
  void bfs(const int root){
    int cnt = 0;
    queue<int> q;
    q.push(root);
    while(!q.empty()){
      int h = q.front(); q.pop();
      for(int i=h; i!=-1; i = heavy_child[i]){
        vid[i] = cnt++;
        vid2idx[vid[i]] = i;
        chain_head[i] = h;
        for(int to : tree[i]) if(to!=par[i] && to!=heavy_child[i]) q.push(to);
      }
    }
  }
public:
  vector<int> tree[SZ];
  int par[SZ], dep[SZ], sub_sz[SZ], heavy_child[SZ], vid[SZ], vid2idx[SZ], chain_head[SZ];
  HLDecomp(const int _n) : n(_n) {
    assert(n<=SZ);
    fill(sub_sz, sub_sz+n, 1);
    fill(heavy_child, heavy_child+n, -1);
  }
  void add_edge(const int u, const int v){
    tree[u].push_back(v);
    tree[v].push_back(u);
  }
  void build(const int root = 0){
    dfs(root);
    bfs(root);
  }

  // 閉区間注意 SegTreeの区間 r+1 必要かも
  void for_each_vertex(int u, int v, const function<void(int,int)> &f) const {
    while(true){
      if(vid[u] > vid[v]) swap(u,v);
      f(max(vid[chain_head[v]], vid[u]), vid[v]);
      if(chain_head[u] != chain_head[v]) v = par[chain_head[v]];
      else break;
    }
  }
  // 辺は子側の頂点で管理
  void for_each_edge(int u, int v, const function<void(int,int)> &f) const {
    while(true){
      if(vid[u] > vid[v]) swap(u,v);
      if(chain_head[u] != chain_head[v]){
        f(vid[chain_head[v]], vid[v]);
        v = par[chain_head[v]];
      }
      else {
        if(u!=v) f(vid[u]+1, vid[v]);
        break;
      }
    }
  }
  int lca(int u, int v) const {
    while(true){
      if(vid[u] > vid[v]) swap(u,v);
      if(chain_head[u] == chain_head[v]) return u;
      v = par[chain_head[v]];
    }
  }
  int distance(const int u, const int v) const {
    return dep[u] + dep[v] - 2*dep[lca(u,v)];
  }
};
