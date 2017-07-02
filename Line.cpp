struct Point{double x,y,mine,high;};
bool comp(Point a,Point b) {return (a.x==b.x ? a.y<b.y : a.x<b.x);}
double dist(Point a, Point b) {return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));}

// 線分p1-p2, p3-p4の交点を求め,resに格納. if文内のu,vの等号で端点を含むかどうか分ける.
bool Intersection(Point p1,Point p2,Point p3,Point p4,Point *res){
  double d = (p2.x-p1.x)*(p4.y-p3.y)-(p2.y-p1.y)*(p4.x-p3.x);
  double u = ((p3.x-p1.x)*(p4.y-p3.y)-(p3.y-p1.y)*(p4.x-p3.x))/d;
  double v = ((p3.x-p1.x)*(p2.y-p1.y)-(p3.y-p1.y)*(p2.x-p1.x))/d;
  if (d == 0 || u<0.0 || u>1.0 || v<0.0 || v>1.0)return false;
  res->x = p1.x + u * (p2.x - p1.x);
  res->y = p1.y + u * (p2.y - p1.y);
  res->high = (p3.mine == p3.high);
  return true;
}
