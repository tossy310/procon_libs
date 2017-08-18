// vector<int> tree[100005];
// の形の木・森に対して同型判定のときに使うコード．
// 根つき木のときはcenterのかわりに根を設定する．

map<vector<int>,int> conv;
int pre[300005];
int val[300005];
bool checked[300005];

// vec上でvを含む木の中心
vector<int> center(int v, vector<int> *vec){
  auto bfs = [&](int d){
    queue<pair<int,int>> q;
    q.push(mp(d,-1));
    val[d] = 0;
    int last = d;
    while(q.size()){
      auto p = q.front(); q.pop();
      int node = p.first;
      int from = p.second;
      checked[node]=true;
      for(auto to : vec[node]) if(to!=from){
        q.push(mp(to,node));
        val[to] = val[node]+1;
        pre[to] = node;
      }
      last = node;
    }
    return last;
  };
  int p = bfs(v);
  int q = bfs(p);

  int d = q;
  rep(i,val[q]/2) d = pre[d];

  if(val[q]%2==1) return {d, pre[d]};
  else return {d};
}

// vec上でiを含む木の同型判定用ナンバリング
pair<int,int> tree2idx(int i, vector<int> *vec){
  auto c = center(i, vec);

  function<int(int,int)> dfs = [&](int d, int from){
    vector<int> ord;
    queue<pair<int,int>> q;
    q.push(mp(d,from));
    while(q.size()){
      auto p = q.front(); q.pop();
      int node = p.first;
      int prev = p.second;
      ord.pb(node);
      pre[node] = prev;
      for(auto to : vec[node]) if(to != prev){
        q.push(mp(to, node));
      }
    }
    reverse(all(ord));
    for(auto node : ord){
      vector<int> v;
      for(auto to : vec[node]) if(to!=pre[node]) v.pb(val[to]);
      sort(all(v));
      if(conv.count(v)==0) conv[v] = conv.size();
      val[node] = conv[v];
    }
    return val[d];
  };

  if(c.size()==1){
    return mp(-1, dfs(c[0],-1));
  }
  else if(c.size()==2){
    int a = dfs(c[0], c[1]);
    int b = dfs(c[1], c[0]);
    return mp(min(a,b), max(a,b));
  }
  else assert(false);
}

vector<pair<int,int>> forest2vec(int sz, vector<int> *vec){
  vector<pair<int,int>> ret;
  fill(checked, checked+sz, false);
  rep(i,sz) if(!checked[i]) ret.pb(tree2idx(i, vec));
  sort(all(ret));
  return ret;
}
