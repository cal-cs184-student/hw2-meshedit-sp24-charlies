#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  {
    std::vector<Vector2D> rst;

    for (int i = 0; i < points.size() - 1; i++)
    {
      Vector2D p1 = points[i];
      Vector2D p2 = points[i + 1];

      rst.push_back(t * p1 + (1 - t) * p2);
    }
    return rst;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    std::vector<Vector3D> rst;

    for (int i = 0; i < points.size() - 1; i++)
    {
      Vector3D p1 = points[i];
      Vector3D p2 = points[i + 1];

      rst.push_back(t * p1 + (1 - t) * p2);
    }
    return rst;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    std::vector<Vector3D> rst = points;
    while (rst.size() > 1)
    {
      rst = BezierPatch::evaluateStep(rst, t);
    }
    return rst[0];
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const
  {
    std::vector<Vector3D> rst;
    for (int i = 0; i < controlPoints.size(); i++)
    {
      rst.push_back(evaluate1D(controlPoints[i], u));
    }

    return evaluate1D(rst, v);
  }

  Vector3D Vertex::normal(void) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
    Vector3D rst = Vector3D(0, 0, 0);

    HalfedgeCIter h = this->halfedge();
    do
    {
      std::vector<Vector3D> vectors;

      HalfedgeCIter vs = h;
      for (int i = 0; i < 3; i++)
      {
        vectors.push_back(vs->vertex()->position);
        vs = vs->next();
      }

      Vector3D v1 = vectors[2] - vectors[0];
      Vector3D v2 = vectors[2] - vectors[1];

      Vector3D n = cross(v1, v2) / cross(v1, v2).norm();

      float area = cross(v1, v2).norm() / 2;

      rst = rst + n * area;

      h = h->twin()->next();

    } while (h != this->halfedge());
    return rst / rst.norm();
  }

  EdgeIter HalfedgeMesh::flipEdge(EdgeIter e0)
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.

    // Halfedges
    HalfedgeIter h0 = e0->halfedge();
    HalfedgeIter h1 = h0->next();
    HalfedgeIter h2 = h1->next();
    HalfedgeIter h3 = h0->twin();
    HalfedgeIter h4 = h3->next();
    HalfedgeIter h5 = h4->next();
    HalfedgeIter h6 = h1->twin();
    HalfedgeIter h7 = h2->twin();
    HalfedgeIter h8 = h5->twin();
    HalfedgeIter h9 = h4->twin();

    // Vertices
    VertexIter v0 = h0->vertex();
    VertexIter v1 = h3->vertex();
    VertexIter v2 = h2->vertex();
    VertexIter v3 = h5->vertex();

    // Edges
    EdgeIter e1 = h1->edge();
    EdgeIter e2 = h2->edge();
    EdgeIter e3 = h4->edge();
    EdgeIter e4 = h5->edge();

    // Faces
    FaceIter f0 = h0->face();
    FaceIter f1 = h3->face();

    // Update
    if (f0->isBoundary() || f1->isBoundary())
    {
      return e0;
    }

    // Halfedges
    h0->next() = h1;
    h0->twin() = h3;
    h0->vertex() = v2;
    h0->edge() = e0;
    h0->face() = f0;

    h1->next() = h2;
    h1->twin() = h8;
    h1->vertex() = v3;
    h1->edge() = e4;
    h1->face() = f0;

    h2->next() = h0;
    h2->twin() = h6;
    h2->vertex() = v1;
    h2->edge() = e1;
    h2->face() = f0;

    h3->next() = h4;
    h3->twin() = h0;
    h3->vertex() = v3;
    h3->edge() = e0;
    h3->face() = f1;

    h4->next() = h5;
    h4->twin() = h7;
    h4->vertex() = v2;
    h4->edge() = e2;
    h4->face() = f1;

    h5->next() = h3;
    h5->twin() = h9;
    h5->vertex() = v0;
    h5->edge() = e3;
    h5->face() = f1;

    h6->next() = h6->next();
    h6->twin() = h2;
    h6->vertex() = h6->vertex();
    h6->edge() = h6->edge();
    h6->face() = h6->face();

    h7->next() = h7->next();
    h7->twin() = h4;
    h7->vertex() = h7->vertex();
    h7->edge() = h7->edge();
    h7->face() = h7->face();

    h8->next() = h8->next();
    h8->twin() = h1;
    h8->vertex() = h8->vertex();
    h8->edge() = h8->edge();
    h8->face() = h8->face();

    h9->next() = h9->next();
    h9->twin() = h5;
    h9->vertex() = h9->vertex();
    h9->edge() = h9->edge();
    h9->face() = h9->face();

    // Vertices
    v0->halfedge() = h5;
    v1->halfedge() = h2;
    v2->halfedge() = h0;
    v3->halfedge() = h3;

    // Edges
    e0->halfedge() = h0;
    e1->halfedge() = h2;
    e2->halfedge() = h4;
    e3->halfedge() = h5;
    e4->halfedge() = h1;

    // Faces
    f0->halfedge() = h0;
    f1->halfedge() = h3;

    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge(EdgeIter e0)
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.

    // Halfedges
    HalfedgeIter h0 = e0->halfedge();
    HalfedgeIter h1 = h0->next();
    HalfedgeIter h2 = h1->next();
    HalfedgeIter h3 = h0->twin();
    HalfedgeIter h4 = h3->next();
    HalfedgeIter h5 = h4->next();
    HalfedgeIter h6 = h1->twin();
    HalfedgeIter h7 = h5->twin();
    HalfedgeIter h8 = h4->twin();
    HalfedgeIter h9 = h2->twin();

    // Vertices
    VertexIter v0 = h0->vertex();
    VertexIter v1 = h1->vertex();
    VertexIter v2 = h2->vertex();
    VertexIter v3 = h5->vertex();

    // Edges
    EdgeIter e1 = h1->edge();
    EdgeIter e2 = h2->edge();
    EdgeIter e3 = h4->edge();
    EdgeIter e4 = h5->edge();

    // Faces
    FaceIter f0 = h0->face();
    FaceIter f1 = h3->face();

    // Create new halfedges
    HalfedgeIter h10 = newHalfedge();
    // HalfedgeIter h11 = newHalfedge();
    HalfedgeIter h12 = newHalfedge();
    HalfedgeIter h13 = newHalfedge();
    HalfedgeIter h14 = newHalfedge();
    HalfedgeIter h15 = newHalfedge();
    // HalfedgeIter h16 = newHalfedge();
    HalfedgeIter h17 = newHalfedge();

    // Create new vertices
    VertexIter v4 = newVertex();

    // Create new edges
    EdgeIter e6 = newEdge();
    EdgeIter e7 = newEdge();
    EdgeIter e8 = newEdge();

    // Create new faces
    FaceIter f2 = newFace();
    FaceIter f3 = newFace();
    FaceIter f4 = newFace();
    FaceIter f5 = newFace();

    // Assign all edges
    h0->setNeighbors(h1, h3, v4, e0, f2);
    h1->setNeighbors(h10, h6, v1, e1, f2);
    h2->setNeighbors(h12, h9, v2, e2, f3);
    h3->setNeighbors(h17, h0, v1, e0, f5);
    h4->setNeighbors(h14, h8, v0, e3, f4);
    h5->setNeighbors(h3, h7, v3, e4, f5);
    h6->setNeighbors(h6->next(), h1, h6->vertex(), h6->edge(), h6->face());
    h7->setNeighbors(h7->next(), h5, h7->vertex(), h7->edge(), h7->face());
    h8->setNeighbors(h8->next(), h4, h8->vertex(), h8->edge(), h8->face());
    h9->setNeighbors(h9->next(), h2, h9->vertex(), h9->edge(), h9->face());

    h10->setNeighbors(h0, h13, v2, e6, f2);
    // h11->setNeighbors(h1, h16, v4, e5, f2);
    h12->setNeighbors(h13, h15, v0, e7, f3);
    h13->setNeighbors(h2, h10, v4, e6, f3);
    h14->setNeighbors(h15, h17, v3, e8, f4);
    h15->setNeighbors(h4, h12, v4, e7, f4);
    // h16->setNeighbors(h17, h11, v1, e5, f5);
    h17->setNeighbors(h5, h14, v4, e8, f5);

    // Assign all vertices
    v0->halfedge() = h4;
    v1->halfedge() = h1;
    v2->halfedge() = h2;
    v3->halfedge() = h5;
    v4->halfedge() = h0;

    // Update the position of v4
    v4->position = (v0->position + v1->position) / 2;
    v4->isNew = true;

    // Assign all edges
    // e5->halfedge() = h11;
    e0->halfedge() = h0;
    e6->halfedge() = h10;
    e7->halfedge() = h12;
    e8->halfedge() = h14;

    e0->isSplit = true;
    e7->isSplit = true;

    e6->isNew = true;
    e8->isNew = true;

    // Assign all faces
    f2->halfedge() = h1;
    f3->halfedge() = h2;
    f4->halfedge() = h4;
    f5->halfedge() = h5;

    // Delete any unused.
    // deleteHalfedge(h0);
    // deleteHalfedge(h3);

    // deleteEdge(e0);

    deleteFace(f0);
    deleteFace(f1);

    return v4;
  }

  void MeshResampler::upsample(HalfedgeMesh &mesh)
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.

    // Compute all the positions of updated old vertices.
    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++)
    {
      Size n = v->degree();

      float u;
      if (n == 3)
      {
        u = 3.0f / 16.0f;
      }
      else
      {
        u = 3.0f / (8.0f * n);
      }
      v->isNew = false;

      // Compute the sum of the neighbor vertices
      Vector3D neighbor_sum;
      HalfedgeIter h = v->halfedge();
      do {
        neighbor_sum += h->twin()->vertex()->position;
        h = h->twin()->next();
      } while (h != v->halfedge());

      v->newPosition = (1.0 - n * u) * v->position + u * neighbor_sum;
      // cout << v->position << " " << v->newPosition << endl;
    }


    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.

    // Compute all the positions of new vertices.
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++)
    {
      // Get A and B vertices.
      VertexIter a = e->halfedge()->vertex();
      VertexIter b = e->halfedge()->twin()->vertex();

      // Get C and D vertices
      VertexIter c = e->halfedge()->next()->next()->vertex();
      VertexIter d = e->halfedge()->twin()->next()->next()->vertex();

      // // Mark as old
      // e->isNew = false;
      e->isNew = false;
      e->isSplit = false;

      // cout << a->position << " " << b->position << " " << c->position << " " << d->position << endl;
      // exit(0);

      // Calculate the new position
      e->newPosition = 3.0 * (a->position + b->position) / 8.0 + 1.0  * (c->position + d->position) / 8.0;
    }

    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)

    // Split all the edges
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++)
    {
      if (!e->isNew && !e->isSplit)
      {
        Vector3D newPosition = e->newPosition;
        VertexIter v = mesh.splitEdge(e);
        v->position = newPosition;
      }
    }
  

    // // 4. Flip any new edge that connects an old and new vertex.
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++)
    {
      if (e->isNew && ((e->halfedge()->vertex()->isNew && !e->halfedge()->twin()->vertex()->isNew) || (!e->halfedge()->vertex()->isNew && e->halfedge()->twin()->vertex()->isNew)))
      {
        mesh.flipEdge(e);
      }
    }

    // // 5. Copy the new vertex positions into final Vertex::position.
    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++)
    {
      if (!v->isNew) {
        v->position = v->newPosition;
      }
      
    }
  }
}
