// TODO Geometry.cpp への完全移行

#define EPS 1e-8

// 浮動小数点誤差考慮
inline double add(double a, double b){
  if(abs(a+b) < EPS*(abs(a) + abs(b))) return 0;
  return a+b;
}

struct Point{
  double x,y;
  Point() {}
  Point(double nx, double ny) : x(nx), y(ny) {}
  inline Point operator + (const Point p){ return Point(add(x, p.x), add(y, p.y)); }
  inline Point operator - (const Point p){ return Point(add(x,-p.x), add(y,-p.y)); }
  inline Point operator * (double d){ return Point(x*d, y*d); }
  inline double dot(const Point p){ return add(x * p.x, y*p.y); }  //内積
  inline double det(const Point p){ return add(x * p.y, -y*p.x); } //外積
  inline double norm(){ return x*x + y*y; }
  inline double dist(const Point & p){ return hypot(x-p.x, y-p.y); }
  inline bool operator < (const Point & p) const {
    if(x != p.x) return x < p.x;
    else return y < p.y;
  }
  inline bool operator == (const Point & p) const {
    return (add(x, -p.x)==0) && (add(y, -p.y)==0);
  }
  friend ostream& operator<<(ostream& os, const Point& p) {
    os << "[" << p.x << "," << p.y << "]";
    return os;
  }
};


// Graham Scan
// 与えらえた点集合の凸包を求める
// 左下の点(Point のoperatorに依存)から反時計回りに
// while内の等号を抜くと辺上の点も含む．
vector<Point> GrahamScan(vector<Point> points){
  int n = points.size();
  sort(all(points));
  int k=0;
  vector<Point> qs(2*n);
  for(int i=0; i<n; qs[k++] = points[i++]){
    while(k>1 && (qs[k-1] - qs[k-2]).det(points[i] - qs[k-1]) <= 0) k--;
  }
  for(int i=n-2,t=k; i>=0; qs[k++] = points[i--]){
    while(k>t && (qs[k-1] -qs[k-2]).det(points[i] - qs[k-1]) <= 0) k--;
  }
  qs.resize(k-1);
  return qs;
}


// キャリパー法 最遠点対探索，最遠距離を返す
// 入力は凸包となっていること．
double caliper(vector<Point> conv){
  int n = conv.size();
  if(n==2) return conv[0].dist(conv[1]);
  int i=0, j=0;
  rep(k, n){
    if(conv[k].x < conv[i].x) i=k;
    if(conv[k].x > conv[j].x) j=k;
  }
  // iが左端，jが右端
  double res=conv[i].dist(conv[j]);
  int si=i, sj=j;
  while(si != j || sj != i){
    if( (conv[(i+1)%n]-conv[i]).det(conv[(j+1)%n]-conv[j]) < 0 ) i = (i+1)%n;
    else j = (j+1)%n;
    res = max(res, conv[i].dist(conv[j]));
  }
  return res;
}

// 点pの直線a(a1-a2)への射影点を返す
Point proj(Point a1, Point a2, Point p) {
  Point a = a2-a1;
  return a1 + a * ( a.dot(p-a1)/a.norm() );
}

// 線分a1-a2上にpがあるか
bool on_line(Point a1, Point a2, Point p){
  Point q1 = a1-p, q2 = a2-p;
  return (q1.det(q2) == 0) && (q1.dot(q2) <= 0);
}

// 線分p1-p2 と 線分p3-p4の交点を求め，resに格納．
// if文内のu,vの等号で端点を含むかどうか分ける．（返り値に影響）
bool line_intersection(Point &p1, Point &p2, Point &p3, Point &p4, Point &res){
  double d =  (p2.x-p1.x)*(p4.y-p3.y) - (p2.y-p1.y)*(p4.x-p3.x);
  double u = ((p3.x-p1.x)*(p4.y-p3.y) - (p3.y-p1.y)*(p4.x-p3.x))/d;
  double v = ((p3.x-p1.x)*(p2.y-p1.y) - (p3.y-p1.y)*(p2.x-p1.x))/d;
  if (u<0.0 || u>1.0 || v<0.0 || v>1.0) return false;
  if(d==0){ // 並行のとき，trueは返せるが交点は不定
    return on_line(p1,p2,p3) || on_line(p1,p2,p4) || on_line(p3,p4,p1) || on_line(p3,p4,p2);
  }
  res.x = p1.x + u * (p2.x - p1.x);
  res.y = p1.y + u * (p2.y - p1.y);
  return true;
}



Point points[50000];
Point conv[50000];
// Graham Scan for POJ, 引数でもとの点の個数，返り値で凸法頂点数
// 左下の点(Point のoperatorに依存)から反時計回りに
// while内の等号を抜くと辺上の点も含む．
int GrahamScan(int n){
  int k=0;
  for(int i=0; i<n; conv[k++] = points[i++]){
    while(k>1 && (conv[k-1] - conv[k-2]).det(points[i] - conv[k-1]) <= 0) k--;
  }
  for(int i=n-2,t=k; i>=0; conv[k++] = points[i--]){
    while(k>t && (conv[k-1] -conv[k-2]).det(points[i] - conv[k-1]) <= 0) k--;
  }
  return k-1;
}
