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
  if(cross(b,c) < -EPS) return -1;  // ccw
  if(dot(b,c)   < -EPS) return +2;  // c--a--b on line
  if(norm(b) < norm(c)) return -2;  // a--b--c on line or a==b
  return 0;                          // a--c--b on line or a==c or b==c
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

double distLP(Point a1, Point a2, Point p){
  return abs(proj(a1, a2, p) - p);
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
