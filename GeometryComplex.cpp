using Point = complex<long>;
const double EPS = 1e-8;
#define X real()
#define Y imag()
namespace std {
  bool operator<(const Point a, const Point b) {
    return a.X != b.X ? a.X < b.X : a.Y < b.Y;
  }
}

// 内積　dot(a,b) = |a||b|cosθ
double dot(Point a, Point b){
  return a.X*b.X + a.Y*b.Y;
}
// 外積　cross(a,b) = |a||b|sinθ
double cross(Point a, Point b){
  return a.X*b.Y - a.Y*b.X;
}


// Graham Scan
// O(N logN)
vector<Point> GrahamScan(vector<Point> points){
  int n = points.size();
  sort(all(points));
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
    // hull.size() >= 2
    auto l = begin(hull);
    auto r = --end(hull);
    if(p.X < l->X || p.X > r->X) return false;

    auto q1 = hull.upper_bound((Point){p.X, INF});
    if( q1==end(hull) ){
      return (--q1)->Y >= p.Y;
    }
    auto q2 = q1--; //dbg(*q1,*q2, q2==end(hull));
    // q1.x <= p.x < q2.x
    return cross(p - *q1, *q2 - *q1) >= -EPS;
  }

  bool inPolygon(Point p){
    bool bu = isCovered(upperHull, p);
    bool bl = isCovered(lowerHull, (Point){p.X, -p.Y});
    return bu==true && bl==true;
  }

  // only for upper hull
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
