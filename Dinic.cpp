// Dinic Maximum Flow algorithm O(E V^2)
template<typename T>
class MaximumFlow{
private:
  struct edge{int to; T cap; int rev;};
  vector<vector<edge> > Graph;
  vector<int> level, iter; //sからの距離,どこまで調べたか
  void bfs(int s){
    fill(begin(level), end(level), -1);
    queue<int> q;
    level[s]=0;
    q.push(s);
    while(!q.empty()){
      int v=q.front(); q.pop();
      for(auto e : Graph[v]){
        if(e.cap>0 && level[e.to]<0){
          level[e.to] = level[v]+1;
          q.push(e.to);
        }
      }
    }
  }
  T dfs(int v, int t, T f){
    if(v==t) return f;
    for(int &i=iter[v]; i<(int)Graph[v].size(); i++){
      auto &e = Graph[v][i];
      if(e.cap>0 && level[v]<level[e.to]){
        T d = dfs(e.to, t, min(f, e.cap));
        if(d>0){
          e.cap -= d;
          Graph[e.to][e.rev].cap += d;
          return d;
        }
      }
    }
    return 0;
  }
public:
  MaximumFlow(int n){
    Graph = vector<vector<edge> >(n, vector<edge>());
    level = vector<int>(n);
    iter = vector<int>(n);
  }
  T max_flow(int s, int t){
    T flow=0;
    while(true){
      bfs(s);
      if(level[t] < 0) break;
      fill(begin(iter), end(iter), 0);
      T f;
      while((f=dfs(s,t,INF)) > 0){
        flow += f;
      }
    }
    return flow;
  }
  void add_edge(int from, int to, T cap){
    int tos = Graph[to].size(), froms = Graph[from].size();
    Graph[from].push_back(((edge){to, cap, tos}));
    Graph[to].push_back(((edge){from, 0, froms}));
  }
}; // END class MaximumFlow



// for POJ
// Dinic Maximum Flow algorithm O(E V^2)
template<typename T>
class MaximumFlow{
private:
  struct edge{int to; T cap; int rev;};
  vector<vector<edge> > Graph;
  vector<int> level, iter; //sからの距離,どこまで調べたか
  void bfs(int s){
    fill(all(level), -1);
    queue<int> q;
    level[s]=0;
    q.push(s);
    while(!q.empty()){
      int v=q.front(); q.pop();
      rep(i, Graph[v].size()){
        edge e = Graph[v][i];
        if(e.cap>0 && level[e.to]<0){
          level[e.to] = level[v]+1;
          q.push(e.to);
        }
      }
    }
  }
  T dfs(int v, int t, T f){
    if(v==t) return f;
    for(int &i=iter[v]; i<(int)Graph[v].size(); i++){
      edge &e = Graph[v][i];
      if(e.cap>0 && level[v]<level[e.to]){
        T d = dfs(e.to, t, min(f, e.cap));
        if(d>0){
          e.cap -= d;
          Graph[e.to][e.rev].cap += d;
          return d;
        }
      }
    }
    return 0;
  }
public:
  MaximumFlow(int n){
    Graph = vector<vector<edge> >(n, vector<edge>());
    level = vector<int>(n);
    iter = vector<int>(n);
  }
  T max_flow(int s, int t){
    T flow=0;
    while(true){
      bfs(s);
      if(level[t] < 0) break;
      fill(all(iter), 0);
      T f;
      while( (f=dfs(s,t,INF)) > 0){
        flow += f;
      }
    }
    return flow;
  }
  void add_edge(int from, int to, T cap){
    int tos = Graph[to].size(), froms = Graph[from].size();
    Graph[from].pb(((edge){to, cap, tos}));
    Graph[to].pb(((edge){from, 0, froms}));
  }
}; // END class MaximumFlow
