// Segment Tree (range min, point update), INF注意
template<typename T>
class SegTree {
public:
  int n;
  vector<T> data;
  SegTree(int n_){
    n=1;
    while(n<n_) n*=2;
    data.resize(2*n, INF);
  }
  inline T query(int l, int r){
    T ret = INF;
    for(l+=n, r+=n; l<r; l/=2, r/=2){
      if(l&1) ret = min(ret, data[l++]);
      if(r&1) ret = min(ret, data[--r]);
    }
    return ret;
  }
  inline void update(int k, T a){
    k+=n;
    data[k]=a;
    while(k>0){
      k = k/2;
      data[k] = min(data[k*2], data[k*2+1]);
    }
  }
  inline T operator[](int idx){ return data[idx+n]; }
};


// RMQ Position (range min position, point update), INF注意
template<typename T>
class SegTree {
public:
  int n;
  vector<T> data;
  vector<int> pos;
  SegTree(int n_){
    n=1;
    while(n<n_) n*=2;
    data.resize(n+1, INF);
    pos.resize(2*n, n);
  }
  inline int query(int l, int r){ // position of min
    int ret = n;
    for(l+=n, r+=n; l<r; l/=2, r/=2){
      if(l&1){
        if(data[ret] >= data[pos[l]]) ret = pos[l];
        l++;
      }
      if(r&1){
        r--;
        if(data[ret] >= data[pos[r]]) ret = pos[r];
      }
    }
    return ret;
  }
  inline void update(int k, T a){
    data[k]=a;
    pos[k+n]=k;
    k+=n;
    while(k>0){
      k = k/2;
      if(data[pos[k*2]] <= data[pos[k*2+1]]){
        pos[k] = pos[k*2];
      }
      else {
        pos[k] = pos[k*2+1];
      }
    }
  }
  inline T operator[](int idx){ return data[idx]; }
};


// Segment Tree (range add, point get), INF注意
template<typename T>
class SegTree {
private:
  int n;
  vector<T> data;
public:
  SegTree(int n_){
    n=1;
    while(n<n_) n*=2;
    data = vector<T>(2*n-1, 0);
  }
  void update(int a, int b, T v, int k=0, int l=0, int r=-1){
    if(r==-1) r=n;
    if(r<=a || b<=l) return;
    if(a<=l && r<=b){
      data[k]+=v;
      return;
    }
    update(a,b,v,k*2+1,l,(l+r)/2);
    update(a,b,v,k*2+2,(l+r)/2,r);
    return;
  }
  T get(int idx){
    idx += n-1;
    T ret = data[idx];
    while( idx>0 ){
      idx = (idx-1) / 2;
      ret += data[idx];
    }
    return ret;
  }
  inline T operator[](int idx){ return data[idx+n-1]; }
};


// Starry Sky Stree (Range Add, Range Max)
// 初期値注意
template<typename T>
class SegTree {
public:
  int n;
  vector<T> segMax, segAdd;
  vector<int> segIdx;
  void _add(int a, int b, T x, int k, int l, int r){
    if(r<=a || b<=l) return;
    if(a<=l && r<=b){segAdd[k]+=x; return;}
    int cl = k*2+1, cr = k*2+2;
    _add(a,b,x,cl,l,(l+r)/2);
    _add(a,b,x,cr,(l+r)/2,r);
    T vl = segMax[cl]+segAdd[cl];
    T vr = segMax[cr]+segAdd[cr];
    if(vl > vr){
      segMax[k] = vl;
      segIdx[k] = segIdx[cl];
    } else {
      segMax[k] = vr;
      segIdx[k] = segIdx[cr];
    }
  }
  pair<T, int> _max(int a, int b, int k, int l, int r){ // <val, idx> (tree上でのidx)
    if(r<=a || b<=l) return mp(-INF, -100);
    if(a<=l && r<=b) return mp(segMax[k]+segAdd[k], segIdx[k]);
    pair<T,int> vl = _max(a,b,k*2+1,l,(l+r)/2);
    pair<T,int> vr = _max(a,b,k*2+2,(l+r)/2,r);
    if(vl.fi > vr.fi) return mp(vl.fi+segAdd[k], vl.se);
    else return mp(vr.fi+segAdd[k], vr.se);
  }
  SegTree(int n_){
    n=1;
    while(n<n_) n*=2;
    segMax.resize(2*n-1);
    fill(all(segMax), 0);
    segAdd.resize(2*n-1);
    fill(all(segAdd), 0);
    segIdx.resize(2*n-1);
    repl(i,n-1,2*n-1) segIdx[i] = i;
    for(int i=n-2; i>=0; i--) segIdx[i] = segIdx[2*i+1];
  }
  inline void add(int a, int b, T x){ _add(a,b,x,0,0,n);} // add x in [a,b)
  inline pair<T,int> getMax(int a,int b){return _max(a,b,0,0,n);} // <max-val, idx> idx はst.n-1を引いたほうがいいかも
};


// Segment Tree (Range Add, Range Sum)
template<typename T>
class SegTree {
private:
  int n;
  vector<T> segAll, segPart;
  void _add(int a, int b, T x, int k, int l, int r){
    if(r<=a || b<=l) return;
    if(a<=l && r<=b){segAll[k]+=x; return;}
    int cl = k*2+1, cr = k*2+2;
    _add(a,b,x,cl,l,(l+r)/2);
    _add(a,b,x,cr,(l+r)/2,r);
    segPart[k] += (min(b,r)-max(a,l))*x;
  }
  T _sum(int a, int b, int k, int l, int r){
    if(r<=a || b<=l) return 0;
    if(a<=l && r<=b) return ((r-l)*segAll[k] + segPart[k]);
    T vl = _sum(a,b,k*2+1,l,(l+r)/2);
    T vr = _sum(a,b,k*2+2,(l+r)/2,r);
    return vl + vr + (min(b,r)-max(a,l))*segAll[k];
  }
public:
  SegTree(int n_){
    n=1;
    while(n<n_) n*=2;
    segAll.resize(2*n-1);
    fill(all(segAll), 0);
    segPart.resize(2*n-1);
    fill(all(segPart), 0);
  }
  inline void add(int a, int b, T x){ _add(a,b,x,0,0,n);} //add x in [a,b)
  inline T getSum(int a,int b){return _sum(a,b,0,0,n);} //sum in [a,b)
};


// Lazy Propagation Segment Tree (Range update, Range sum)
class SegTree {
public:
  int n;
  vector<long> data, sum; // data: 区間内の値，sum:区間内の実際の和
  inline void lazy_eval(int k, int l, int r){
    if(data[k]>0){ // 最新の更新が下に伝搬していないとき
      sum[k] = data[k]*(r-l);
      if(k<n-1){
        data[2*k+1]=data[k];
        data[2*k+2]=data[k];
      }
      data[k]=0;
    }
  }
  void _update(int a, int b, long v, int k, int l, int r){
    lazy_eval(k,l,r);
    if(r<=a || b<=l) return;
    if(a<=l && r<=b){
      data[k] = v;
      lazy_eval(k,l,r);
      return;
    }
    int m=(l+r)/2;
    _update(a,b,v,k*2+1,l,m);
    _update(a,b,v,k*2+2,m,r);
    sum[k] = sum[k*2+1] + sum[k*2+2];
    return;
  }
  // [a,b) をvに更新
  void update(int a, int b, long v){ _update(a,b,v,0,0,n); }

  long _sum(int a, int b, int k, int l, int r){
    lazy_eval(k,l,r);
    if(r<=a || b<=l) return 0;
    if(a<=l && r<=b) return sum[k];
    int m=(l+r)/2;
    long vl = _sum(a,b,k*2+1,l,m);
    long vr = _sum(a,b,k*2+2,m,r);
    sum[k] = sum[k*2+1] + sum[k*2+2];
    return vl+vr;
  }
  // [a,b) の合計
  long query(int a,int b){ return _sum(a,b,0,0,n); }

  SegTree(int n_){
    n=1;
    while(n<n_) n*=2;
    data = vector<long>(2*n-1, 0);
    sum = vector<long>(2*n-1, 0);
  }
};


// ノードに区間を持つセグメント木
// グローバルの配列vecからn個値をとってきている．
class SegTree{
public:
  int n;
  vector<int> data[1<<18];
  SegTree(int n_){
    n=1;
    while(n<n_) n*=2;
    rep(i,n_){
      data[n-1+i].pb(vec[i]);
    }
    for(int i=n-2; i>=0; i--){
      int il = 2*i+1, ir = 2*i+2;
      data[i].resize(data[il].size()+data[ir].size());
      merge(all(data[il]), all(data[ir]), data[i].begin());
    }
  }
  // 区間[a,b)でx以下の個数
  int query(int a, int b, int x, int k, int l, int r){
    if(r<=a || b<=l) return 0;
    if(a<=l && r<=b){
      return upper_bound(all(data[k]), x) - data[k].begin();
    }
    int vl = query(a,b,x,k*2+1,l,(r+l)/2);
    int vr = query(a,b,x,k*2+2,(r+l)/2,r);
    return vl+vr;
  }
  int query(int a, int b, int x){
    return query(a,b,x,0,0,n);
  }
};
