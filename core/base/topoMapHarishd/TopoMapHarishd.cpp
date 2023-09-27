#include <TopoMapHarishd.h>
#include <vector>
#include <set>
#include <queue>
#include <map>
#include <utility>
#include <algorithm>
#include <Geometry.h>
#include <string>
#include <cmath>

ttk::TopoMapHarishd::TopoMapHarishd() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("TopoMapHarishd");
}
ttk::TopoMapHarishd::~TopoMapHarishd() = default;

std::vector<Point> ttk::TopoMapHarishd::project(std::vector<double> &data, int dimension) {
  std::vector<std::pair<int, int> > edges;
  std::vector<double> weights;

  log("computing emst");
  this->emst(data,dimension,edges,weights);
  log("placing points");
  return this->placePoints(edges,weights, data.size()/dimension);
}

void ttk::TopoMapHarishd::emst(std::vector<double> &data, int dim, std::vector<std::pair<int, int> > &edges, std::vector<double> &weights) {
  int npts = data.size() / dim;
  for(int i = 0;i < npts;i ++) {
    for(int j = 0;j < i; j ++) {
      double dx = data[3*i]-data[3*j];
      double dy = data[3*i+1]-data[3*j+1];
      double dz = data[3*i+2]-data[3*j+2];
      double dist = sqrt(dx*dx+dy*dy+dz*dz);
      weights.push_back(dist);
      edges.push_back({i,j});
    }
  }
  return;
}

std::vector<Point> ttk::TopoMapHarishd::placePoints(const std::vector<std::pair<int, int> > &edges, const std::vector<double> &weights, int n) {
  if(edges.size() != weights.size()) {
    std::cerr << "edges and lengths don't match" << std::endl;
    exit(0);
  }
  comps = DisjointSets<int>(n);
  compMap.clear();
  compMap.resize(n);
  verts.resize(n);

  for(size_t i = 0;i < compMap.size();i ++) {
    verts[i].p = Point(0,0);
    verts[i].id = i;
    compMap[i].vertices.push_back(i);
    compMap[i].hull.push_back(verts[i].p);
    compMap[i].hull.push_back(verts[i].p);
  }
  std::vector<int> order = sortEdges(edges,weights);

  for(size_t _i = 0;_i < order.size();_i ++) {
    size_t i = order[_i];
    int p1 = edges[i].first;
    int p2 = edges[i].second;

    int c1 = comps.find(p1);
    int c2 = comps.find(p2);

    if(c1 == c2) {
      continue;
      //std::cerr << "Error!!! MST edge belongs to the same component!!!" << std::endl;
      //exit(0);
    }

    Component &comp1 = compMap[c1];
    Component &comp2 = compMap[c2];

    Component comp = mergeComponents(comp1, comp2, p1, p2, weights[i]);
    comps.merge(c1, c2);
    int c = comps.find(c1);
    compMap[c] = comp;
  }
  std::vector<Point> pts;
  for(size_t i = 0;i < verts.size();i ++) {
    pts.push_back(verts[i].p);
  }
  return pts;
}

Transformation ttk::TopoMapHarishd::alignHull(const Polygon &hull, const Point &p, bool topEdge) {
  int v = -1;
  double d2;
  for(size_t i = 0;i < hull.size()-1;i ++) {
    double d = distance2(hull[i],p);
    if(v == -1) {
      d2 = d;
      v = i;
    } else {
      if(d2 > d) {
        d2 = d;
        v = i;
      }
    }
  }

  // hull is ordered cloclwise by boost
  // v is the required vertex
  Point v1, v2;
  if(topEdge) {
    // make v,v+1 the top edge of the hull, s.t. v = (0,0);
    v1 = hull[v];
    v2 = hull[v+1];
  } else {
    // make v,v-1 the bottom edge of the hull, s.t. v = (0,0);
    if(v == 0) {
      v = hull.size() - 1;
    }
    v1 = hull[v];
    v2 = hull[v-1];
  }

  Transformation trans;
  // translate
  trans.tx = -hull[v].x;
  trans.ty = -hull[v].y;
  // rotate
  if(hull.size() > 2) {
    findAngle(v1,v2,trans);
  } else {
    trans.sin = 0;
    trans.cos = 1;
  }
  return trans;
}

inline Point transform(const Point &p, const Transformation &t, double yOffset) {
  double x = p.x + t.tx;
  double y = p.y + t.ty;

  double xx = x * t.cos - y * t.sin;
  double yy = x * t.sin + y * t.cos;

  yy += yOffset;
  return Point(xx,yy);
}

void ttk::TopoMapHarishd::transformComponent(const Component &c, const Transformation &t, double yOffset) {
  for(size_t i = 0;i < c.vertices.size();i ++) {
    int vin = c.vertices[i];
    this->verts[vin].p = transform(this->verts[vin].p,t, yOffset);
  }
}

Component ttk::TopoMapHarishd::mergeComponents(Component &c1, Component &c2, int v1, int v2, double length) {
  // compute the merged component
  Component merged;
  merged.vertices.clear();
  merged.vertices.insert(merged.vertices.end(),c1.vertices.begin(),c1.vertices.end());
  merged.vertices.insert(merged.vertices.end(),c2.vertices.begin(),c2.vertices.end());

  merged.hull.clear();
  if(length > 0) {
    // transform the vertices of the two components appropriately
    Transformation t1 = alignHull(c1.hull,verts[v1].p,true); // aligned w.r.t. top edge. so stays in the bottom
    this->transformComponent(c1,t1,0);

    Transformation t2 = alignHull(c2.hull,verts[v2].p,false); // aligned w.r.t. bottom edge. so offset should be added
    this->transformComponent(c2,t2,length);

    // compute the hull of the merged component
    Polygon pts;
    for(size_t i = 0;i < c1.hull.size() - 1;i ++) {
      pts.push_back(transform(c1.hull[i],t1,0));
    }
    for(size_t i = 0;i < c2.hull.size() - 1;i ++) {
      pts.push_back(transform(c2.hull[i],t2,length));
    }
    computeConvexHull(pts,merged.hull);
  } else {
    // support for points sharing the same coordinates (MST edge length = 0)
    if(c1.hull.size() != 2 || c2.hull.size() != 2) {
      std::cerr << "Error!!! hull cannot have more than one point when edge lenght is 0!!!" << std::endl;
      exit(0);
    }
    merged.hull.push_back(c2.hull[0]);
    merged.hull.push_back(c2.hull[1]);
  }

  return merged;
}

void ttk::TopoMapHarishd::log(std::string str) {
  if(verbose) {
    std::cout << str << std::endl;
  }
}



