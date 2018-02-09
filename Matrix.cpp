typedef vector<vector<long>> mat;

// return A*B
mat mat_mul(const mat &A, const mat &B){
  int n=A.size(), m=B[0].size(), l=B.size();
  mat ret(n, vector<long>(m, 0));
  rep(i,n) rep(k,l) if(A[i][k]!=0) rep(j,m){
    (ret[i][j] += A[i][k] * B[k][j]) %= MOD;
  }
  return ret;
}

// A^p
mat mat_pow(const mat &A, long p){
  int n = A.size();
  mat tmp(A), ret(n, vector<long>(n,0));
  rep(i,n) ret[i][i] = 1;
  while(p>0){
    if(p&1) ret = mat_mul(tmp, ret);
    tmp = mat_mul(tmp, tmp);
    p /= 2;
  }
  return ret;
}
// Aが零行列で零行列を返さないので注意．（0^xの定義にもよる）

#define EPS 1e-9
// solve Ax=b O(N^3)
vector<double> gaussJordan(const vector<vector<double>> &A, const vector<double> &b){
  int n=A.size();
  vector<vector<double>> B(n, vector<double>(n+1));
  rep(i,n) rep(j,n) B[i][j] = A[i][j];
  rep(i,n) B[i][n] = b[i];

  rep(i,n){
    int pivot=i;
    rep(j,i,n) if(abs(B[j][i]) > abs(B[pivot][i])) pivot = j;
    swap(B[i], B[pivot]);

    if(abs(B[i][i]) < EPS) return vector<double>();

    rep(j,i+1,n+1) B[i][j] /= B[i][i];
    rep(j,n) if(i!=j) rep(k,i+1,n+1) B[j][k] -= B[j][i] * B[i][k];
  }
  vector<double> x(n);
  rep(i,n) x[i] = B[i][n];
  return x;
}


// Matrix 2*2 (require ModInt)
struct mat {
  Int a,b,c,d;
  mat() : a(1), b(0), c(0), d(1) {}
  mat(Int a0, Int a1, Int a2, Int a3) : a(a0), b(a1), c(a2), d(a3) {}
  mat &operator *= (const mat &p){
    Int aa = a*p.a + b*p.c;
    Int bb = a*p.b + b*p.d;
    Int cc = c*p.a + d*p.c;
    Int dd = c*p.b + d*p.d;
    a = aa; b = bb; c = cc; d = dd;
    return *this;
  }
  mat operator * (const mat &p) const {
    return mat(*this) *= p;
  }
  mat inverse() const {
    Int id = (a*d - b*c).inverse();
    return mat(id*d, -id*b, -id*c, id*a);
  }
  friend ostream &operator << (ostream &os, const mat &p) {
    os << p.a << " " << p.b << " " << p.c << " " << p.d;
    return os;
  }
};


// for POJ
#define MAT_N 101
long mat[MAT_N][MAT_N];

// a * b = dst O(N^3)
long tmp[MAT_N][MAT_N];
void mat_mul(long a[][MAT_N], long b[][MAT_N], long dst[][MAT_N]){
  fill(tmp[0], tmp[MAT_N], 0);
  rep(i,MAT_N) rep(k,MAT_N) if(a[i][k]!=0) rep(j,MAT_N){
    tmp[i][j] += a[i][k] * b[k][j];
  }
  rep(i,MAT_N) rep(j,MAT_N) dst[i][j] = tmp[i][j];
}

// a^n = dst O(M^3 log N)
long tmp_pow[MAT_N][MAT_N];
void mat_pow(long a[][MAT_N], long n, long dst[][MAT_N]){
  rep(i,MAT_N) rep(j,MAT_N) tmp_pow[i][j] = a[i][j];
  fill(dst[0], dst[MAT_N], 0);
  rep(i,MAT_N) dst[i][i]=1;
  while(n>0){
    if(n&1) mat_mul(tmp_pow, dst, dst);
    mat_mul(tmp_pow, tmp_pow, tmp_pow);
    n /= 2;
  }
}
