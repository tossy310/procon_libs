using Point = complex<double>;
const double EPS = 1e-8;
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

// 内積 dot(a,b) = |a||b|cosθ
double dot(Point a, Point b){ return a.X*b.X + a.Y*b.Y; }
// 外積 cross(a,b) = |a||b|sinθ
double cross(Point a, Point b){ return a.X*b.Y - a.Y*b.X; }

// AB からみて AC がどの方向にあるか
int ccw(Point a, Point b, Point c){
  b -= a;  c -= a;
  if(cross(b,c) >  EPS) return +1;  // ccw
  if(cross(b,c) < -EPS) return -1;  // cw
  if(dot(b,c)   < -EPS) return +2;  // c--a--b on line
  if(norm(b) < norm(c)) return -2;  // a--b--c on line or a==b
  return 0;                          // a--c--b on line or a==c or b==c
}

// 交差判定
bool isecLP(Point a1, Point a2, Point b){
  return abs(ccw(a1, a2, b)) != 1;  // return EQ(cross(a2-a1, b-a1), 0);
}
bool isecLL(Point a1, Point a2, Point b1, Point b2){
  return !isecLP(a2-a1, b2-b1, 0) || isecLP(a1, b1, b2);
}
bool isecLS(Point a1, Point a2, Point b1, Point b2) {
  return cross(a2-a1, b1-a1) * cross(a2-a1, b2-a1) < EPS;
}
bool isecSP(Point a1, Point a2, Point b){
  return !ccw(a1, a2, b);
}
bool isecSS(Point a1, Point a2, Point b1, Point b2){
  return ccw(a1, a2, b1)*ccw(a1, a2, b2) <= 0 && ccw(b1, b2, a1)*ccw(b1, b2, a2) <= 0;
}


// 点pの直線a1-a2への射影点
Point proj(Point a1, Point a2, Point p){
  return a1 + dot(a2-a1, p-a1)/norm(a2-a1) * (a2-a1);
}

// 距離
double distLP(Point a1, Point a2, Point p){
  return abs(proj(a1, a2, p) - p);
}
double distLL(Point a1, Point a2, Point b1, Point b2){
  return isecLL(a1, a2, b1, b2) ? 0 : distLP(a1, a2, b1);
}
double distLS(Point a1, Point a2, Point b1, Point b2){
  return isecLS(a1, a2, b1, b2) ? 0 : min(distLP(a1, a2, b1), distLP(a1, a2, b2));
}
double distSP(Point a1, Point a2, Point p){
  Point r = proj(a1, a2, p);
  if(isecSP(a1, a2, r)) return abs(r-p);
  return min(abs(a1-p), abs(a2-p));
}
double distSS(Point a1, Point a2, Point b1, Point b2){
  if(isecSS(a1, a2, b1, b2)) return 0;
  return min(min(distSP(a1, a2, b1), distSP(a1, a2, b2)), min(distSP(b1, b2, a1), distSP(b1, b2, a2)));
}
double distLC(Point a1, Point a2, Point c, double r){
  return max(distLP(a1, a2, c) - r, 0.0);
}
double distSC(Point a1, Point a2, Point c, double r){ // not verified
  double dSqr1 = norm(c-a1), dSqr2 = norm(c-a2);
  bool b1 = dSqr1 < r*r, b2 = dSqr2 < r*r;
  if(b1 ^ b2) return 0; // 交差
  if(b1 & b2) return r - sqrt(max(dSqr1, dSqr2)); // 内包. 場合により0
  return max(distSP(a1, a2, c) - r, 0.0);
}
double distCC(Point a, double ar, Point b, double br){
  double d = abs(a-b);
  return GE(d, abs(ar-br)) ? max(d-ar-br, 0.0) : abs(ar-br) - d;
}

// 交点
Point crosspointLL(Point a1, Point a2, Point b1, Point b2){
  double d1 = cross(b2-b1, b1-a1);
  double d2 = cross(b2-b1, a2-a1);
  if(EQ(d1, 0) && EQ(d2, 0)) return a1;  // same line
  assert(!EQ(d2, 0));
  return a1 + d1/d2 * (a2-a1);
}
vector<Point> crosspointLC(Point a1, Point a2, Point c, double r){ // not verified
  vector<Point> ps;
  Point ft = proj(a1, a2, c);
  if(!GE(r*r, norm(ft-c))) return ps;
  Point dir = sqrt(max(r*r - norm(ft-c), 0.0)) / abs(a2-a1) * (a2-a1);
  ps.pb(ft + dir);
  if(!EQ(r*r, norm(ft-c))) ps.pb(ft - dir);
  return ps;
}
vector<Point> crosspointCC(Point a, double ar, Point b, double br){ // not verified
  vector<Point> ps;
  Point ab = b-a;
  double d = abs(ab);
  double crL = (norm(ab) + ar*ar - br*br) / (2*d);
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
int inConvex(Point p, const vector<Point>& conv){
  int n = conv.size();
  int dir = ccw(conv[0], conv[1], p);
  rep(i,n){
    int ccwc = ccw(conv[i], conv[(i+1)%n], p);
    if (!ccwc) return 2;  // 境界上
    if (ccwc != dir) return 0;
  }
  return 1;
}

// 凸多角形の内部判定 O(logn)
// 領域内部なら1、境界上なら2、外部なら0
int inConvex(Point p, const vector<Point>& conv){
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
  double cr = cross(conv[l] - p, conv[r] - p);
  return EQ(cr, 0) ? 2 : cr < 0 ? 0 : 1;
}


// キャリパー法 最遠点対探索，最遠距離を返す
// 入力は凸包となっていること．
double caliper(vector<Point> conv){
  int n = conv.size();
  if(n==2) return abs(conv[0] - conv[1]);
  int i=0, j=0;
  rep(k, n){
    if(conv[k].x < conv[i].x) i=k;
    if(conv[k].x > conv[j].x) j=k;
  }
  // iが左端，jが右端
  double res = abs(conv[i] - conv[j]);
  int si=i, sj=j;
  while(si != j || sj != i){
    if( cross(conv[(i+1)%n]-conv[i], conv[(j+1)%n]-conv[j]) < 0 ) i = (i+1)%n;
    else j = (j+1)%n;
    res = max(res, abs(conv[i]-conv[j]));
  }
  return res;
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

// 点pから円aへの接線の接点
VP tangentPoints(P a, D ar, P p) {
  VP ps;
  D sin = ar / abs(p-a);
  if (!LE(sin, 1)) return ps;  // ここでNaNも弾かれる
  D t = M_PI_2 - asin(min(sin, 1.0));
  ps.push_back(                 a + (p-a)*polar(sin, t));
  if (!EQ(sin, 1)) ps.push_back(a + (p-a)*polar(sin, -t));
  return ps;
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

// 点aと点bを通り、半径がrの円の中心を返す
VP circlesPointsRadius(P a, P b, D r) {
  VP cs;
  P abH = (b-a)*0.5;
  D d = abs(abH);
  if (d == 0 || d > r) return cs;  // 必要なら !LE(d,r) として円1つになる側へ丸める
  D dN = sqrt(r*r - d*d);          // 必要なら max(r*r - d*d, 0) とする
  P n = abH * P(0,1) * (dN / d);
  cs.push_back(a + abH + n);
  if (dN > 0) cs.push_back(a + abH - n);
  return cs;
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

// 多角形の符号付面積
D area(const VP& ps) {
  D a = 0;
  rep (i, ps.size()) a += cross(ps[i], ps[(i+1) % ps.size()]);
  return a / 2;
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
