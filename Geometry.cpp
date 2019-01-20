using Double = double;
using Point = complex<Double>;
const Double EPS = 1e-8;
#define X real()
#define Y imag()
#define LE(n,m) ((n) < (m) + EPS)
#define GE(n,m) ((n) + EPS > (m))
#define EQ(n,m) (abs((n)-(m)) < EPS)
namespace std {
  bool operator<(const Point a, const Point b){
    return a.X != b.X ? a.X < b.X : a.Y < b.Y;
  }
}

// =========== In case only integer arithmetics ===========
using Point = pair<long long, long long>;
#define X first
#define Y second
namespace std {
  Point operator+(const Point &a, const Point &b){
    return {a.X + b.X, a.Y + b.Y};
  }
  Point operator-(const Point &a, const Point &b){
    return {a.X - b.X, a.Y - b.Y};
  }
  bool operator<(const Point &a, const Point &b){
    return a.X != b.X ? a.X < b.X : a.Y < b.Y;
  }
}
// ===========
/*
  p = polar(1, r) : abs(p) = 1, arg(p) = r [rad] を生成。オイラーの公式より、これを掛ければ原点中心に r 回転。
*/

// 内積 dot(a,b) = |a||b|cosθ
Double dot(const Point &a, const Point &b){ return a.X*b.X + a.Y*b.Y; }
// 外積 cross(a,b) = |a||b|sinθ
Double cross(const Point &a, const Point &b){ return a.X*b.Y - a.Y*b.X; }

// AB からみて AC がどの方向にあるか
int ccw(const Point &a, const Point &b, const Point &c){
  const Point db = b - a;
  const Point dc = c - a;
  if(cross(db,dc) >  EPS) return +1;  // ccw
  if(cross(db,dc) < -EPS) return -1;  // cw
  if(dot(db,dc)   < -EPS) return +2;  // c--a--b on line
  if(norm(db) < norm(dc)) return -2;  // a--b--c on line or a==b
  return 0;                           // a--c--b on line or a==c or b==c
}

// 交差判定
bool isecLP(const Point &a1, const Point &a2, const Point &b){
  return abs(ccw(a1, a2, b)) != 1;  // return EQ(cross(a2-a1, b-a1), 0);
}
bool isecLL(const Point &a1, const Point &a2, const Point &b1, const Point &b2){
  return !isecLP(a2-a1, b2-b1, 0) || isecLP(a1, b1, b2);
}
bool isecLS(const Point &a1, const Point &a2, const Point %b1, const Point &b2){
  return cross(a2-a1, b1-a1) * cross(a2-a1, b2-a1) < EPS;
}
bool isecSP(const Point &a1, const Point &a2, const Point &b){
  return !ccw(a1, a2, b);
}
bool isecSS(const Point &a1, const Point &a2, const Point &b1, const Point &b2){
  return ccw(a1, a2, b1)*ccw(a1, a2, b2) <= 0 && ccw(b1, b2, a1)*ccw(b1, b2, a2) <= 0;
}


// 点pの直線a1-a2への射影点
Point proj(const Point &a1, const Point &a2, const Point &p){
  return a1 + dot(a2-a1, p-a1)/norm(a2-a1) * (a2-a1);
}

// 距離
Double distLP(const Point &a1, const Point &a2, const Point &p){
  return abs(proj(a1, a2, p) - p);
}
Double distLL(const Point &a1, const Point &a2, const Point &b1, const Point &b2){
  return isecLL(a1, a2, b1, b2) ? 0 : distLP(a1, a2, b1);
}
Double distLS(const Point &a1, const Point &a2, const Point &b1, const Point &b2){
  return isecLS(a1, a2, b1, b2) ? 0 : min(distLP(a1, a2, b1), distLP(a1, a2, b2));
}
Double distSP(const Point &a1, const Point &a2, const Point &p){
  Point r = proj(a1, a2, p);
  if(isecSP(a1, a2, r)) return abs(r-p);
  return min(abs(a1-p), abs(a2-p));
}
Double distSS(const Point &a1, const Point &a2, const Point &b1, const Point &b2){
  if(isecSS(a1, a2, b1, b2)) return 0;
  return min(min(distSP(a1, a2, b1), distSP(a1, a2, b2)), min(distSP(b1, b2, a1), distSP(b1, b2, a2)));
}
Double distLC(const Point &a1, const Point &a2, const Point &c, const Double r){
  return max(distLP(a1, a2, c) - r, 0.0);
}
Double distSC(const Point &a1, const Point &a2, const Point &c, const Double r){ // not verified
  Double dSqr1 = norm(c-a1), dSqr2 = norm(c-a2);
  bool b1 = dSqr1 < r*r, b2 = dSqr2 < r*r;
  if(b1 ^ b2) return 0; // 交差
  if(b1 & b2) return r - sqrt(max(dSqr1, dSqr2)); // 内包. 場合により0
  return max(distSP(a1, a2, c) - r, 0.0);
}
Double distCC(const Point &a, const Double ar, const Point &b, const Double br){
  Double d = abs(a-b);
  return GE(d, abs(ar-br)) ? max(d-ar-br, 0.0) : abs(ar-br) - d;
}

// 交点
Point crosspointLL(const Point &a1, const Point &a2, const Point &b1, const Point &b2){
  Double d1 = cross(b2-b1, b1-a1);
  Double d2 = cross(b2-b1, a2-a1);
  if(EQ(d1, 0) && EQ(d2, 0)) return a1;  // same line
  assert(!EQ(d2, 0));
  return a1 + d1/d2 * (a2-a1);
}
vector<Point> crosspointLC(const Point &a1, const Point &a2, const Point &c, const Double r){
  vector<Point> ps;
  Point ft = proj(a1, a2, c);
  if(!GE(r*r, norm(ft-c))) return ps;
  Point dir = sqrt(max(r*r - norm(ft-c), 0.0)) / abs(a2-a1) * (a2-a1);
  ps.pb(ft + dir);
  if(!EQ(r*r, norm(ft-c))) ps.pb(ft - dir);
  return ps;
}
vector<Point> crosspointCC(const Point &a, const Double ar, const Point &b, const Double br){ // not verified
  vector<Point> ps;
  Point ab = b-a;
  Double d = abs(ab);
  Double crL = (norm(ab) + ar*ar - br*br) / (2*d);
  if(EQ(d, 0) || ar < abs(crL)) return ps;
  Point abN = ab * Point(0, sqrt(ar*ar - crL*crL) / d);
  Point cp = a + crL/d * ab;
  ps.pb(cp + abN);
  if(!EQ(norm(abN), 0)) ps.pb(cp - abN);
  return ps;
}


// Graham Scan
// O(N logN)
vector<Point> GrahamScan(vector<Point> points){
  uniq(points);
  int n=points.size();
  int k=0;
  vector<Point> qs(2*n);
  for(int i=0; i<n; qs[k++] = points[i++]){
    while(k>1 && cross(qs[k-1] - qs[k-2], points[i] - qs[k-1]) <= 0) k--;
  }
  for(int i=n-2,t=k; i>=0; qs[k++] = points[i--]){
    while(k>t && cross(qs[k-1] -qs[k-2], points[i] - qs[k-1]) <= 0) k--;
  }
  qs.resize(k-1);
  return qs;
}


// 凸包の内部判定 O(n)
// 領域内部なら1、境界上なら2、外部なら0
int inConvex(const Point &p, const vector<Point>& conv){
  int n = conv.size();
  int dir = ccw(conv[0], conv[1], p);
  rep(i,n){
    int ccwc = ccw(conv[i], conv[(i+1)%n], p);
    if (!ccwc) return 2;  // 境界上
    if (ccwc != dir) return 0;
  }
  return 1;
}

// 凸包の内部判定 O(logn)
// 領域内部なら1、境界上なら2、外部なら0
int inConvex(const Point &p, const vector<Point>& conv){
  int n = conv.size();
  Point g = (conv[0] + conv[n/3] + conv[n*2/3])/3.0;
  if (g == p) return 1;
  Point gp = p - g;

  int l=0, r=n;
  while(r-l>1){
    int mid = (l+r)/2;
    Point gl = conv[l] - g;
    Point gm = conv[mid] - g;
    if(cross(gl, gm) > 0){
      if(cross(gl, gp) >= 0 && cross(gm, gp) <= 0) r = mid;
      else l = mid;
    }
    else {
      if(cross(gl, gp) <= 0 && cross(gm, gp) >= 0) l = mid;
      else r = mid;
    }
  }
  r %= n;
  Double cr = cross(conv[l] - p, conv[r] - p);
  return EQ(cr, 0) ? 2 : cr < 0 ? 0 : 1;
}

// 多角形の符号付面積
Double area(const vector<Point>& ps) {
  Double a = 0;
  rep(i, ps.size()){
    a += cross(ps[i], ps[(i+1)%ps.size()]);
  }
  return a/2;
}


// キャリパー法 最遠点対探索，最遠距離を返す
// 入力は凸包となっていること．
Double caliper(const vector<Point> &conv){
  int n = conv.size();
  if(n==2) return abs(conv[0] - conv[1]);
  int i=0, j=0;
  rep(k, n){
    if(conv[k].x < conv[i].x) i=k;
    if(conv[k].x > conv[j].x) j=k;
  }
  // iが左端，jが右端
  Double res = abs(conv[i] - conv[j]);
  int si=i, sj=j;
  while(si != j || sj != i){
    if( cross(conv[(i+1)%n]-conv[i], conv[(j+1)%n]-conv[j]) < 0 ) i = (i+1)%n;
    else j = (j+1)%n;
    res = max(res, abs(conv[i]-conv[j]));
  }
  return res;
}


// 点aと点bを通り、半径がrの円の中心を返す
vector<Point> circlesPointsRadius(const Point &a, const Point &b, const Ddouble r){
  vector<Point> cs;
  Ppoint abH = (b-a)*0.5;
  Double d = abs(abH);
  if(d == 0 || d > r) return cs;  // a==b  or  abs(a-b) > 2r
  Double dN = sqrt(r*r - d*d);
  Point n = abH * P(0,1) * (dN / d);
  cs.push_back(a + abH + n);
  if(dN > 0) cs.push_back(a + abH - n);
  return cs;
}

// 点pから円aへの接線の接点
vector<Point> tangentPoints(const Point &a, const Double ar, const Point &p) {
  vector<Point> ps;
  Double sin = ar / abs(p-a);
  if (!LE(sin, 1)) return ps;  // ここでNaNも弾かれる
  Double t = M_PI_2 - asin(min(sin, 1.0));
  ps.push_back(a + (p-a)*polar(sin, t));
  if (!EQ(sin, 1)) ps.push_back(a + (p-a)*polar(sin, -t));
  return ps;
}

// ConvexHull with Vertex addition
// O((N+Q)logN)
class AddableConvexHull {
public:
  set<Point> upperHull, lowerHull;
  AddableConvexHull(){}

  bool isCovered(set<Point> &hull, Point p){
    if(hull.size()==0) return false;
    if(hull.size()==1){
      Point q = *begin(hull);
      return q.X == p.X && q.Y >= p.Y;
    }
    auto l = begin(hull);
    auto r = --end(hull);
    if(p.X < l->X || p.X > r->X) return false;

    auto q1 = hull.upper_bound((Point){p.X, INF});
    if( q1==end(hull) ){
      return (--q1)->Y >= p.Y;
    }
    auto q2 = q1--;
    // q1.x <= p.x < q2.x
    return cross(p - *q1, *q2 - *q1) >= -EPS;
  }

  bool inPolygon(Point p){
    bool bu = isCovered(upperHull, p);
    bool bl = isCovered(lowerHull, (Point){p.X, -p.Y});
    return bu==true && bl==true;
  }

  void add(set<Point> &hull, Point p){
    if(hull.size()==0){
      hull.insert(p);
      return;
    }
    if(hull.size()==1){
      if(begin(hull)->X == p.X){
        Point q = *begin(hull);
        hull.clear();
        hull.insert(max(p,q));
      }
      else {
        hull.insert(p);
      }
      return;
    }
    if(isCovered(hull, p)) return;

    // left
    auto q1 = hull.upper_bound((Point){p.X, INF});
    if(q1 != begin(hull)){
      q1--;
      while(hull.size()>1 && q1!=begin(hull)){
        auto q2 = q1--;
        if( cross(*q2 - p, *q1 - p) > -EPS ) break;
        hull.erase(q2);
      }
    }
    // right
    q1 = hull.lower_bound((Point){p.X, -INF});
    if(q1 != end(hull)){
      while(hull.size()>1 && q1 != --end(hull)){
        auto q2 = q1++;
        if( cross(*q1 - p, *q2 - p) > -EPS ) break;
        hull.erase(q2);
      }
    }

    hull.insert(p);
  }

  void add(Point p){
    add(upperHull, p);
    add(lowerHull, (Point){p.X, -p.Y});
  }
};



// 余弦定理
// △ABC において、a = BC, b = CA, c = AB としたとき
// a^2 = b^2 + c^2 ? 2bc cos ∠CAB
//
// ヘロンの公式
// 3辺の長さがa,b,cである三角形の面積T
// T = sqrt{ s(s-a)(s-b)(s-c) }, s = (a+b+c)/2
//
//
//
// 以下は色々なとこからコピペしたものを貼ってあるだけで信憑性低め

// ベクトルpをベクトルbに射影したベクトルを計算する
inline P proj(const P& p, const P& b) {
    return b*inp(p,b)/norm(b);
}
// 点pから直線lに引いた垂線の足となる点を計算する
inline P perf(const L& l, const P& p) {
    L m = {l.pos - p, l.dir};
    return (p + (m.pos - proj(m.pos, m.dir)));
}
// 線分sを直線bに射影した線分を計算する
inline L proj(const L& s, const L& b) {
     return (L){perf(b, s.pos), proj(s.dir, b.dir)};
}

// 2円の共通接線。返される各直線に含まれる頂点は円との接点となる
vector<L> tangentLines(P a, D ar, P b, D br) {
  vector<L> ls;
  D d = abs(b-a);
  rep (i,2) {
    D sin = (ar - (1-i*2)*br) / d;
    if (!LE(sin*sin, 1)) break;
    D cos = sqrt(max(1 - sin*sin, 0.0));
    rep (j,2) {
      P n = (b-a) * P(sin, (1-j*2)*cos) / d;
      ls.push_back(L(a + ar*n, b + (1-i*2)*br*n));
      if (cos < EPS) break;  // 重複する接線を無視（重複していいならこの行不要）
    }
  }
  return ls;
}

// 2円の共通面積
double cc_area(const C& c1, const C& c2) {
    double d = abs(c1.p - c2.p);
    if (c1.r + c2.r <= d + EPS) {
        return 0.0;
    } else if (d <= abs(c1.r - c2.r) + EPS) {
        double r = c1.r <? c2.r;
        return r * r * PI;
    } else {
        double rc = (d*d + c1.r*c1.r - c2.r*c2.r) / (2*d);
        double theta = acos(rc / c1.r);
        double phi = acos((d - rc) / c2.r);
        return c1.r*c1.r*theta + c2.r*c2.r*phi - d*c1.r*sin(theta);
    }
}

// 三角形の外心。点a,b,cは同一線上にあってはならない
P circumcenter(P a, P b, P c) {
  a = (a-c)*0.5;
  b = (b-c)*0.5;
  return c + crosspointLL(a, a*P(1,1), b, b*P(1,1));
}


// 点aと点bを通り、直線lに接する円の中心
VP circlesPointsTangent(P a, P b, P l1, P l2) {
  P n = (l2-l1) * P(0,1);
  P m = (b-a) * P(0,0.5);
  D rC = dot((a+b)*0.5-l1, n);
  D qa = norm(n)*norm(m) - dot(n,m)*dot(n,m);
  D qb = -rC * dot(n,m);
  D qc = norm(n)*norm(m) - rC*rC;
  D qd = qb*qb - qa*qc;  // qa*k^2 + 2*qb*k + qc = 0

  VP cs;
  if (qd < -EPS) return cs;
  if (EQ(qa, 0)) {
    if (!EQ(qb, 0)) cs.push_back((a+b)*0.5 - m * (qc/qb/2));
    return cs;
  }
  D t = -qb/qa;
  cs.push_back(              (a+b)*0.5 + m * (t + sqrt(max(qd, 0.0))/qa));
  if (qd > EPS) cs.push_back((a+b)*0.5 + m * (t - sqrt(max(qd, 0.0))/qa));
  return cs;
}

// 凸多角形クリッピング
// たぶんpsは反時計回り，直線で切り取られた左側がreturnされていそう．たぶん．
VP convexCut(const VP& ps, P a1, P a2) {
  int n = ps.size();
  VP ret;
  rep(i,n) {
    int ccwc = ccw(a1, a2, ps[i]);
    if (ccwc != -1) ret.push_back(ps[i]);
    int ccwn = ccw(a1, a2, ps[(i + 1) % n]);
    if (ccwc * ccwn == -1) ret.push_back(crosspointLL(a1, a2, ps[i], ps[(i + 1) % n]));
  }
  return ret;
}

// 多角形の幾何学的重心
P centroid(const VP& ps) {
  int n = ps.size();
  D aSum = 0;
  P c;
  rep (i, n) {
    D a = cross(ps[i], ps[(i+1) % n]);
    aSum += a;
    c += (ps[i] + ps[(i+1) % n]) * a;
  }
  return 1 / aSum / 3 * c;
}
