template<typename T>
class SparseTable {
public:
  vector<int> lg_tbl;
  vector<T> tbl[20];
  int n;
  T (*op)(const T&, const T&);
  SparseTable(const vector<T> &vec, T (*_op)(const T&, const T&)) : op(_op){
    n = vec.size();
    lg_tbl.resize(n+1, 0);
    int v = 0;
    rep(i,3,n+1){
      if((1<<v)*2<i) v++;
      lg_tbl[i] = v;
    }
    tbl[0] = vec;
    for(int k=1; (1<<k)<=n; k++){
      tbl[k].resize(n);
      rep(i,n){
        tbl[k][i] = op(tbl[k-1][i], tbl[k-1][min(n-1,i+(1<<(k-1)))]);
      }
    }
  }
  T query(int l, int r){  // [l,r)
    int lg = lg_tbl[r-l];
    return op(tbl[lg][l], tbl[lg][r-(1<<lg)]);
  }
};

// Disjoint version.
template<typename T>
class SparseTable {
public:
  vector<T> tbl[20];
  vector<T> tbr[20];
  int n;
  T (*op)(const T&, const T&);
  SparseTable(const vector<T> &v, T (*_op)(const T&, const T&)) : op(_op){
    // build table O(NlogN)
    n = v.size();
    tbl[0].resize(n+1);
    tbr[0].resize(n+1);
    rep(i,n) tbl[0][i] = tbr[0][i] = v[i];
    for(int k=1; (1<<k)<=n; k++){
      tbl[k].resize(n+1); tbr[k].resize(n+1);
      int mask = (1<<k)-1;
      rep(i,n){
        if((i&mask) == 0) tbr[k][i] = v[i];
        else tbr[k][i] = op(tbr[k][i-1], v[i]);
      }
      for(int i=n-1; i>=0; i--){
        if(((i+1)&mask) == 0) tbl[k][i] = v[i];
        else tbl[k][i] = op(tbl[k][i+1], v[i]);
      }
    }
  }
  T query(int l, int r) const { // [l,r)
    if(l+1 == r) return tbl[0][l];
    r--; // [l,r]
    int k = 31 - __builtin_clz(l^r); // 2^k <= l^r < 2^{k+1}
    return op(tbl[k][l], tbr[k][r]);
  }
};
