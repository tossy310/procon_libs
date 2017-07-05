/*
  3
1 0 4
  2
  5
*/
int dice[6][6] ={ // [top][front] -> right
  {-1, 2, 4, 1, 3,-1},
  { 3,-1, 0, 5,-1, 2},
  { 1, 5,-1,-1, 0, 4},
  { 4, 0,-1,-1, 5, 1},
  { 2,-1, 5, 0,-1, 3},
  {-1, 3, 1, 4, 2,-1}
};


/*
  2
1 0 4
  3
  5
*/
int dice[6][6] ={ // [top][front] -> right
  {-1, 3, 1, 4, 2,-1},
  { 2,-1, 5, 0,-1, 3},
  { 4, 0,-1,-1, 5, 1},
  { 1, 5,-1,-1, 0, 4},
  { 3,-1, 0, 5,-1, 2},
  {-1, 2, 4, 1, 3,-1}
};


/* 以下の座標系で動かすときのdx,dy,df
 ┼─→ y
 │
 ↓
 x

[k]  1
    ↑
3 ←  → 2
    ↓
    0
*/
const int dx[] = {1,-1,0,0}, dy[]={0,0,1,-1};
function<pair<int,int>(int,int)> df[] = { // return <top,front>
  [&](int t, int f){return mp(5-f, t);},
  [&](int t, int f){return mp(f, 5-t);},
  [&](int t, int f){return mp(5-dice[t][f], f);},
  [&](int t, int f){return mp(dice[t][f], f);}
};
