/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/accel.h>
#include <Eigen/Geometry>
#include <queue>
#include <chrono>

template<typename Base, typename T>
inline bool instanceof(const T *ptr) {
  return dynamic_cast<const Base*>(ptr) != nullptr;
}
NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

// TODO: NEED TO TEST THIS FUNCTION
Accel::Node* Accel::buildTree(BoundingBox3f boundingBox, std::vector<uint32_t> triangles,int treeDepth) {
  if(triangles.empty()) {
    return Singleton;
  }

  if(triangles.size() <= 10 || treeDepth > 20) {
    BlackNode* temp = new BlackNode {};
    temp->triangleIndices = triangles;
    return temp;
  }

  Point3f center = boundingBox.getCenter();
  Point3f max = boundingBox.max;
  Point3f min = boundingBox.min;

  BoundingBox3f subBoxes[8] = {
   BoundingBox3f(min,center),
   BoundingBox3f(Point3f(center(0,0),min(1,0),min(2,0)),Point3f(max(0,0),center(1,0),center(2,0))),
   BoundingBox3f(Point3f(min(0,0),center(1,0),min(2,0)),Point3f(center(0,0),max(1,0),center(2,0))),
   BoundingBox3f(Point3f(center(0,0),center(1,0),min(2,0)),Point3f(max(0,0),max(1,0),center(2,0))),
   BoundingBox3f(Point3f(min(0,0),min(1,0),center(2,0)),Point3f(center(0,0),center(1,0),max(2,0))),
   BoundingBox3f(Point3f(center(0,0),min(1,0),center(2,0)),Point3f(max(0,0),center(1,0),max(2,0))),
   BoundingBox3f(Point3f(min(0,0),center(1,0),center(2,0)),Point3f(center(0,0),max(1,0),max(2,0))),
   BoundingBox3f(center,max)
  };

  std::vector<uint32_t> childTriangles[8] = {};
  MatrixXu m_F = m_mesh->getIndices();
  MatrixXf m_V = m_mesh->getVertexPositions();

  for(uint32_t idx = 0; idx < triangles.size(); idx++) { //For all triangles for this subnode

    uint32_t index = triangles[idx];
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);
    BoundingBox3f triangleBox = BoundingBox3f();
    triangleBox.expandBy(p0);
    triangleBox.expandBy(p1);
    triangleBox.expandBy(p2);

    for(uint32_t jdx = 0; jdx < 8; jdx++) { // For all cubes
      if(triangleBox.overlaps(subBoxes[jdx])) {
        childTriangles[jdx].push_back(index);
      }
    }
  }

  GrayNode* parent = new GrayNode {};

  for(uint32_t idx = 0; idx < 8; idx++) {
    parent->children[idx] = buildTree(subBoxes[idx],childTriangles[idx],treeDepth++);
  }
  return parent;
}

void Accel::build() {
  
  cout << std::to_string(sizeof(BlackNode)) <<"\n";
  cout << std::to_string(sizeof(WhiteNode))<<"\n";
  cout << std::to_string(sizeof(GrayNode))<<"\n";
  cout << std::to_string(sizeof(Node))<<"\n";
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  std::vector<uint32_t> triangleIndices = std::vector<uint32_t>();
  for (uint32_t i=0; i< m_mesh->getTriangleCount();i++) {
    triangleIndices.push_back(i);
  }
  m_tree = buildTree(m_bbox,triangleIndices,0);

  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> time_span = t2 - t1;
  std::cout << "It took me " << time_span.count() << " milliseconds.\n";
  std::queue<Node *> traversal = std::queue<Node *>();
  traversal.push(m_tree);
  uint32_t grayNodeCount = 0;
  uint32_t blackNodeCount = 0;
  uint32_t averageTriangles = 0;
  while(!traversal.empty()){
    Node* currNode = traversal.front();
    traversal.pop();
    if(instanceof<GrayNode>(currNode)) { // GrayNode
      GrayNode* grayNode = (GrayNode *) currNode;
      grayNodeCount++;
      for(uint32_t i=0;i<8;i++) {
        traversal.push(grayNode->children[i]);
      }
    }else if(instanceof<BlackNode>(currNode)) { // GrayNode
      BlackNode* blackNode = (BlackNode *) currNode;
      blackNodeCount++;
      averageTriangles += blackNode->triangleIndices.size();
    }
  }
  cout << "There are " << grayNodeCount << "number of Interior nodes\n";
  cout << "There are " << blackNodeCount << "number of Leaf nodes\n";
  cout << "The average number of triangles is  " <<  averageTriangles / blackNodeCount << "\n";

}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    std::queue<std::pair<Node *,BoundingBox3f>> traversal = std::queue<std::pair<Node *,BoundingBox3f>>();
    traversal.push(std::make_pair(m_tree,m_mesh->getBoundingBox()));
    while(!traversal.empty()){
      std::pair<Node *,BoundingBox3f> pair = traversal.front();
      traversal.pop();

      Node* currNode = pair.first;
      BoundingBox3f boundingBox = pair.second;

      if(boundingBox.rayIntersect(ray_)) {
        if(instanceof<GrayNode>(currNode)) { // GrayNode
          GrayNode* grayNode = (GrayNode *) currNode;
          Point3f center = boundingBox.getCenter();
          Point3f max = boundingBox.max;
          Point3f min = boundingBox.min;
          BoundingBox3f subBoxes[8] = {
           BoundingBox3f(min,center),
           BoundingBox3f(Point3f(center(0,0),min(1,0),min(2,0)),Point3f(max(0,0),center(1,0),center(2,0))),
           BoundingBox3f(Point3f(min(0,0),center(1,0),min(2,0)),Point3f(center(0,0),max(1,0),center(2,0))),
           BoundingBox3f(Point3f(center(0,0),center(1,0),min(2,0)),Point3f(max(0,0),max(1,0),center(2,0))),
           BoundingBox3f(Point3f(min(0,0),min(1,0),center(2,0)),Point3f(center(0,0),center(1,0),max(2,0))),
           BoundingBox3f(Point3f(center(0,0),min(1,0),center(2,0)),Point3f(max(0,0),center(1,0),max(2,0))),
           BoundingBox3f(Point3f(min(0,0),center(1,0),center(2,0)),Point3f(center(0,0),max(1,0),max(2,0))),
           BoundingBox3f(center,max)
          };
          for(uint32_t i=0; i < 8; i++) {
            if(subBoxes[i].rayIntersect(ray)) {
              traversal.push(std::make_pair(grayNode->children[i],subBoxes[i]));
            }
          }

        }else if(instanceof<BlackNode>(currNode)) {
          BlackNode* blackNode = (BlackNode *) currNode;
          for (uint32_t idx = 0; idx < blackNode->triangleIndices.size(); ++idx) {
            uint32_t index = blackNode->triangleIndices[idx];
            float u, v, t;
            if (m_mesh->rayIntersect(index, ray, u, v, t)) {
              if (shadowRay)
                return true;
              ray.maxt = its.t = t;
              its.uv = Point2f(u, v);
              its.mesh = m_mesh;
              f = index;
              foundIntersection = true;
            }
    }
        }else if(instanceof<WhiteNode>(currNode)){
        }else{
        }
      }
   }

    /* Brute force search through all triangles
    for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
        float u, v, t;
        if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
            An intersection was found! Can terminate
               immediately if this is a shadow ray query
            if (shadowRay)
                return true;
            ray.maxt = its.t = t;
            its.uv = Point2f(u, v);
            its.mesh = m_mesh;
            f = idx;
            foundIntersection = true;
        }
    }
    */

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END

