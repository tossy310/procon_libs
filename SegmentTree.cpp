template<typename T>
class SegTree {
public:
  int n;
  T e;
  T (*op)(const T&, const T&);
  vector<T> data;
  SegTree(int m, T _e, T (*_op)(const T&, const T&)) : e(_e), op(_op){
    n=1;
    while(n<m) n*=2;
    data.resize(2*n, e);
  }
  SegTree(const vector<T> &v, T _e, T (*_op)(const T&, const T&)) : e(_e), op(_op){
    n=1;
    while(n<(int)v.size()) n*=2;
    data.resize(2*n, e);
    rep(i,v.size()) data[i+n] = v[i];
    for(int i=n-1; i>0; i--) data[i] = op(data[i*2], data[i*2+1]);
  }
  T query(int l, int r) const {
    T vl = e, vr = e;
    for(l+=n, r+=n; l<r; l/=2, r/=2){
      if(l&1) vl = op(vl, data[l++]);
      if(r&1) vr = op(data[--r], vr);
    }
    return op(vl,vr);
  }
  void update(int k, T a){
    k+=n;
    data[k]=a;
    while(k>0){
      k = k/2;
      data[k] = op(data[k*2], data[k*2+1]);
    }
  }
  inline T operator[](int idx) const { return data[idx+n]; }
};


// Segment Tree (range min, point update), INF注意
template<typename T>
class SegTree {
public:
  int n;
  vector<T> data;
  SegTree(int m){
    n=1;
    while(n<m) n*=2;
    data.resize(2*n, INF);
  }
  T query(int l, int r){
    T ret = INF;
    for(l+=n, r+=n; l<r; l/=2, r/=2){
      if(l&1) ret = min(ret, data[l++]);
      if(r&1) ret = min(ret, data[--r]);
    }
    return ret;
  }
  void update(int k, T a){
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
  SegTree(int m){
    n=1;
    while(n<m) n*=2;
    data.resize(n+1, INF);
    pos.resize(2*n, n);
  }
  int query(int l, int r){ // position of min
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
  void update(int k, T a){
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


// RAQ Segment Tree (range add, point get)
template<typename T>
class SegTree {
public:
  int n;
  vector<T> data;
  SegTree(int m){
    n=1;
    while(n<m) n*=2;
    data = vector<T>(2*n, 0);
  }
  void update(int l, int r, T v){
    for(l+=n, r+=n; l<r; l/=2, r/=2){
      if(l&1) data[l++] += v;
      if(r&1) data[--r] += v;
    }
  }
  T get(int k){
    k += n;
    T ret = data[k];
    while(k>0){
      k = k/2;
      ret += data[k];
    }
    return ret;
  }
};


// Starry Sky Stree (Range Add, Range Max) with index
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
    rep(i,n-1,2*n-1) segIdx[i] = i;
    for(int i=n-2; i>=0; i--) segIdx[i] = segIdx[2*i+1];
  }
  SegTree(const vector<T> &v){
    int n_ = v.size();
    n=1;
    while(n<n_) n*=2;
    segMax.resize(2*n-1);
    segAdd.resize(2*n-1, 0);
    segIdx.resize(2*n-1);
    rep(i,n_) segMax[n+i-1] = v[i];
    for(int i=n-2; i>=0; i--) segMax[i] = max(segMax[2*i+1], segMax[2*i+2]);
    rep(i,n-1,2*n-1) segIdx[i] = i;
    for(int i=n-2; i>=0; i--){
      if(segMax[2*i+1] > segMax[2*i+2]) segIdx[i] = 2*i+1;
      else segIdx[i] = 2*i+2;
    }
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
// TODO we assume updated values should be non-negative, thus -1 is a sentinel.
template<typename T>
class SegTree {
public:
  int n;
  vector<T> lazy, val;
  // lazy: uniform value for the range (not propageted), val: actual total sum value of the range
  inline void lazy_eval(int k, int l, int r){
    if(lazy[k]>=0){
      val[k] = lazy[k]*(r-l);
      if(k<n){
        lazy[2*k] = lazy[k];
        lazy[2*k+1] = lazy[k];
      }
      lazy[k] = -1;
    }
  }
  void update(int a, int b, T v, int k, int l, int r){
    lazy_eval(k,l,r);
    if(r<=a || b<=l) return;
    if(a<=l && r<=b){
      lazy[k] = v;
      lazy_eval(k,l,r);
      return;
    }
    int m = (l+r)/2;
    update(a,b,v,k*2,l,m);
    update(a,b,v,k*2+1,m,r);
    val[k] = val[k*2] + val[k*2+1];
    return;
  }
  inline void update(int a, int b, T v){ update(a,b,v,1,0,n); }
  T sum(int a, int b, int k, int l, int r){
    lazy_eval(k,l,r);
    if(r<=a || b<=l) return 0;
    if(a<=l && r<=b) return val[k];
    int m = (l+r)/2;
    T vl = sum(a,b,k*2,l,m);
    T vr = sum(a,b,k*2+1,m,r);
    val[k] = val[k*2] + val[k*2+1];
    return vl + vr;
  }
  inline T sum(int a,int b){ return sum(a,b,1,0,n); }
  SegTree(int n_){
    n=1;
    while(n<n_) n*=2;
    lazy = vector<T>(2*n, -1);
    val = vector<T>(2*n, 0);
  }
  void init(vector<T> &v){
    rep(i,v.size()) val[i+n] = v[i];
    for(int i=n-1; i>0; i--) val[i] = val[i*2] + val[i*2+1];
  }
};


// 点add,区間sum,区間累積和min
template<typename T>
class SegTree {
public:
  int n;
  vector<T> data, accMin;
  SegTree(){}
  SegTree(int m){
    n=1;
    while(n<m) n*=2;
    data.resize(2*n, 0);
    accMin.resize(2*n, INF);
  }
  T getSum(int l, int r){
    T ret = 0;
    for(l+=n, r+=n; l<r; l/=2, r/=2){
      if(l&1) ret += data[l++];
      if(r&1) ret += data[--r];
    }
    return ret;
  }
  T getAccMin(int l, int r){ // <sum, accMin>
    pair<T,T> vl = mp(0,INF), vr = mp(0,INF);
    for(l+=n, r+=n; l<r; l/=2, r/=2){
      if(l&1){
        vl = mp(vl.fi + data[l], min(vl.se, vl.fi+accMin[l]));
        l++;
      }
      if(r&1){
        r--;
        vr = mp(data[r] + vr.fi, min(accMin[r], data[r] + vr.se));
      }
    }
    return min(vl.se, vl.fi + vr.se);
  }
  void update(int k, T a){
    k+=n;
    data[k]=a;
    accMin[k]=a;
    while(k>0){
      k = k/2;
      data[k] = data[k*2] + data[k*2+1];
      accMin[k] = min(accMin[k*2], data[k*2] + accMin[k*2+1]);
    }
  }
  inline T operator[](int idx){ return data[idx+n]; }
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
