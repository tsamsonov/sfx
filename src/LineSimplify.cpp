#include <Rcpp.h>
#include <fstream>
using namespace Rcpp;

double get_distance(double x1, double y1, double x2, double y2) {
  return sqrt(pow(x1-x2, 2) + pow(y1-y2, 2));
}

double get_area(double x1, double y1, double x2, double y2, double x3, double y3){
  double a = get_distance(x1, y1, x2, y2);
  double b = get_distance(x2, y2, x3, y3);
  double c = get_distance(x3, y3, x1, y1);
  
  double p = 0.5 * (a + b + c);
  
  // fabs is important for degenerate triangles, where one (p - ...) could
  // be a small number below zero due to the rounding error
  return sqrt(fabs(p * (p-a) * (p-b) * (p-c)));
}

// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool onSegment(double x1, double y1, double x2, double y2, double x3, double y3)
{
  if (x2 <= std::max(x1, x3) && x2 >= std::min(x1, x3) &&
      y2 <= std::max(y1, y3) && y2 >= std::min(y1, y3))
    return true;
  
  return false;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are collinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(double x1, double y1, double x2, double y2, double x3, double y3)
{
  int val = (y2 - y1) * (x3 - x2) -
    (x2 - x1) * (y3 - y2);
  
  if (val == 0) return 0;  // colinear
  
  return (val > 0)? 1: 2; // clock or counterclock wise
}

// Check intersection 
bool doIntersect(double x11, double y11, double x12, double y12, double x21, double y21, double x22, double y22)
{
  
  int o1 = orientation(x11, y11, x21, y21, x12, y12);
  int o2 = orientation(x11, y11, x22, y22, x12, y12);
  int o3 = orientation(x21, y21, x12, y12, x22, y22);
  int o4 = orientation(x21, y21, x11, y11, x22, y22);
  
  // General case
  if ((o1 != o2) && (o3 != o4)){
    return true;
  }
  
  if (o1 == 0 && onSegment(x11, y11, x21, y21, x12, y12)) return true;
  if (o2 == 0 && onSegment(x11, y11, x22, y22, x12, y12)) return true;
  if (o3 == 0 && onSegment(x21, y21, x12, y12, x22, y22)) return true;
  if (o4 == 0 && onSegment(x21, y21, x11, y11, x22, y22)) return true;
  
  return false; // Doesn't fall in any of the above cases
}

// Check whether the point could be deleted 
bool can_erase(NumericVector x, NumericVector y, int k){
  int npts = x.size();
  
  double x1 = x[k-1];
  double y1 = y[k-1];
  double x2 = x[k+1];
  double y2 = y[k+1];
  
  for (int i = 0; i < k-3; i++){
    bool intersects = doIntersect(x1, y1, x2, y2, x[i], y[i], x[i+1], y[i+1]);
    if(intersects){
      return false;
    }
  }
  for (int i = k+2; i < npts-1; i++){
    bool intersects = doIntersect(x1, y1, x2, y2, x[i], y[i], x[i+1], y[i+1]);
    if(intersects){
      return false;
    }
  }
  
  return true;
}

// [[Rcpp::export]]
NumericVector get_areas(NumericMatrix coords) {
  
  int nareas = coords.nrow() - 2;
  
  if (nareas < 1) {
    return NumericVector::create(0);
  }
  
  NumericVector areas(nareas);
  
  for(int i = 1; i <= nareas; i++) {
    areas[i-1] = get_area(coords(i-1, 0), coords(i-1, 1), 
                          coords(i, 0),   coords(i, 1), 
                          coords(i+1, 0), coords(i+1, 1));
  }
  
  return areas;
}

// [[Rcpp::export]]
int get_ruler(NumericMatrix coords, double dist){
  NumericVector x = coords( _ , 0 );
  NumericVector y = coords( _ , 1 );
  int n = x.size();
  int i = 0;
  int j = 1;
  int k = 0;
  while ((i < n) && (j < n)){
    if((get_distance(x[i], y[i], x[j], y[j]) >= dist)){
      i = j;
      k++;
    }
    j++;
  }
  return k;
}

// [[Rcpp::export]]
NumericMatrix insert_points(NumericMatrix coords, NumericMatrix add){
  
  NumericVector x = coords( _ , 0 );
  NumericVector y = coords( _ , 1 );
  
  NumericVector ax = add( _ , 0 );
  NumericVector ay = add( _ , 1 );
  
  int n = x.size();
  int k = ax.size();
  
  for (int i = 0; i < n-1; i++){
    for (int j = 0; j < k; j++){
      if(onSegment(x[i], y[i], ax[j], ay[j], x[i+1], y[i+1])){
        x.insert(i+1, ax[j]);
        y.insert(i+1, ay[j]);
        ax.erase(j);
        ay.erase(j);
        i++;
        n++;
        k--;
        break;
      }
    }
  }
  
  NumericMatrix output(n, 2);
  output(_,0) = x;
  output(_,1) = y;
  
  return output;
}

// [[Rcpp::export]]
NumericMatrix contract_edges(NumericMatrix coords, double tol){
  
  NumericVector x = coords( _ , 0 );
  NumericVector y = coords( _ , 1 );
  
  int n = x.size();
  
  for (int i = 0; i < n-1; i++){
    if(get_distance(x[i], y[i], x[i+1], y[i+1]) <= tol){
      
      if(i > 0){
        if(i < (n-2)){
          x[i] = 0.5 * (x[i] + x[i+1]);
          y[i] = 0.5 * (y[i] + y[i+1]);
          x.erase(i+1);
          y.erase(i+1);
        } else {
          x.erase(i);
          y.erase(i);
        }
      } else {
        x.erase(i+1);
        y.erase(i+1);
      }
      
      --i;
      --n;
    }
  }
  
  NumericMatrix output(n, 2);
  output(_,0) = x;
  output(_,1) = y;
  
  return output;
}

// [[Rcpp::export]]
int get_nbends(NumericMatrix coords){
  int nbends = 0;
  
  NumericVector x = coords( _ , 0 );
  NumericVector y = coords( _ , 1 );
  
  return nbends;
}

// [[Rcpp::export]]
NumericMatrix simplify_vw(NumericMatrix coords, double minarea){
  
  NumericVector areas = get_areas(coords);
  int nareas = areas.size();
  
  NumericVector x = coords( _ , 0 );
  NumericVector y = coords( _ , 1 );
  
  int index = which_min(areas);
  double area = areas[index];
  
  while (area <= minarea){
    
    if(can_erase(x, y, index + 1)){
      
      areas.erase(index);
      x.erase(index + 1);
      y.erase(index + 1);
      
      nareas--;
      
      if(nareas > 0){
        
        if(index < nareas) {
          areas[index] = get_area(x[index], y[index], 
                                  x[index + 1], y[index + 1], 
                                  x[index + 2], y[index + 2]);        
        }
        
        if(index > 0) {
          areas[index-1] = get_area(x[index - 1], y[index - 1],
                                    x[index], y[index],
                                    x[index + 1], y[index + 1]);
        }
        
      } else break;
    } else {
      areas[index] = minarea + 1;
    }

    index = which_min(areas);
    area = areas[index];
    
  }
  
  NumericMatrix output(nareas + 2, 2);
  output(_,0) = x;
  output(_,1) = y;
  
  return output;
  
}