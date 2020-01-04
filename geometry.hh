#pragma once

#define NO_SURFACE_FIT

#include <array>
#include <iostream>
#include <list>
#include <memory>
#include <vector>

namespace Geometry {

const double epsilon = 1.0e-8;

class Vector2D;
class Vector3D;
class BSCurve;

using Point2D = Vector2D;
using Point3D = Vector3D;
using DoubleVector = std::vector<double>;
using Vector2DVector = std::vector<Vector2D>;
using VectorVector = std::vector<Vector3D>;
using Point2DVector = std::vector<Point2D>;
using PointVector = std::vector<Point3D>;
using CurveVector = std::vector<std::shared_ptr<BSCurve>>;

class Vector2D {
public:
  // Constructors
  Vector2D();
  Vector2D(double x, double y);

  // Assignments
  Vector2D &operator+=(const Vector2D &v);
  Vector2D &operator-=(const Vector2D &v);
  Vector2D &operator*=(double x);
  Vector2D &operator/=(double x);

  // Coordinates
  double *data();
  const double *data() const;
  double &operator[](size_t i);
  const double &operator[](size_t i) const;

  // Arithmetic
  Vector2D operator-() const;
  Vector2D operator+(const Vector2D &v) const;
  Vector2D operator-(const Vector2D &v) const;
  double operator*(const Vector2D &v) const;
  Vector2D operator*(double x) const;
  Vector2D operator/(double x) const;

  // Other
  double norm() const;
  double normSqr() const;
  Vector2D &normalize();

private:
  std::array<double, 2> v_;
};

std::ostream &operator<<(std::ostream &os, const Vector2D &v);
std::istream &operator>>(std::istream &is, Vector2D &v);

class Vector3D {
public:
  // Constructors
  Vector3D();
  Vector3D(double x, double y, double z);

  // Assignments
  Vector3D &operator+=(const Vector3D &v);
  Vector3D &operator-=(const Vector3D &v);
  Vector3D &operator*=(double x);
  Vector3D &operator/=(double x);

  // Coordinates
  double *data();
  const double *data() const;
  double &operator[](size_t i);
  const double &operator[](size_t i) const;

  // Arithmetic
  Vector3D operator-() const;
  Vector3D operator+(const Vector3D &v) const;
  Vector3D operator-(const Vector3D &v) const;
  Vector3D operator^(const Vector3D &v) const;
  double operator*(const Vector3D &v) const;
  Vector3D operator*(double x) const;
  Vector3D operator/(double x) const;

  // Other
  double norm() const;
  double normSqr() const;
  Vector3D &normalize();

private:
  std::array<double, 3> v_;
};

std::ostream &operator<<(std::ostream &os, const Vector3D &v);
std::istream &operator>>(std::istream &is, Vector3D &v);

class Matrix3x3 {
public:
  // Special matrices
  static Matrix3x3 identity();
  static Matrix3x3 rotation(const Vector3D &axis, double angle);

  const double &operator()(size_t i, size_t j) const;
  double &operator()(size_t i, size_t j);

  // Arithmetic
  Matrix3x3 operator+(const Matrix3x3 &m) const;
  Matrix3x3 &operator+=(const Matrix3x3 &m);
  Matrix3x3 operator*(double x) const;
  Matrix3x3 &operator*=(double x);
  Vector3D operator*(const Vector3D &v) const;
  Matrix3x3 operator*(const Matrix3x3 &m) const;
  Matrix3x3 &operator*=(const Matrix3x3 &m);

  Matrix3x3 inverse() const;

private:
  std::array<double, 9> m_;
};

class BSBasis {
public:
  // Constructors
  BSBasis();
  BSBasis(const BSBasis &) = default;
  BSBasis(size_t degree, const DoubleVector &knots);
  BSBasis &operator=(const BSBasis &) = default;

  // Properties
  size_t degree() const;
  void setDegree(size_t degree);
  const DoubleVector &knots() const;
  DoubleVector &knots();

  // Utilities
  void reverse();
  void normalize();
  size_t findSpan(double u) const;
  size_t findSpanWithMultiplicity(double u, size_t &multi) const;
  void basisFunctions(size_t i, double u, DoubleVector &coeff) const;
  void basisFunctionsAll(size_t i, double u, std::vector<DoubleVector> &coeff) const;

private:
  size_t p_;
  DoubleVector knots_;
};

class BSCurve {
public:
  // Constructors
  BSCurve();
  BSCurve(const BSCurve &) = default;
  BSCurve(const PointVector &cpts);
  BSCurve(size_t degree, const DoubleVector &knots, const PointVector &cpts);
  BSCurve &operator=(const BSCurve &) = default;

  // Evaluation
  Point3D eval(double u) const;
  Point3D eval(double u, size_t nr_der, VectorVector &der) const;

  // Coordinates
  const PointVector &controlPoints() const;
  PointVector &controlPoints();

  // Parameterization
  const BSBasis &basis() const;
  void reverse();
  void normalize();

  // Other
  double arcLength(double from, double to) const;
  BSCurve insertKnot(double u, size_t r) const;
  DoubleVector intersectWithPlane(const Point3D &p, const Vector3D &n) const;

private:
  void derivativeControlPoints(size_t d, size_t r1, size_t r2, std::vector<PointVector> &dcp) const;
  BSCurve insertKnot(double u, size_t k, size_t s, size_t r) const;

  size_t n_;
  BSBasis basis_;
  PointVector cp_;
};

class TriMesh {
public:
  using Triangle = std::array<size_t, 3>;

  // Mesh building
  void clear();
  void resizePoints(size_t n);
  void setPoints(const PointVector &pv);
  void addTriangle(size_t a, size_t b, size_t c);
  void setTriangles(const std::list<Triangle> &tl);
  TriMesh &append(const TriMesh &other);

  // I/O
  Point3D &operator[](size_t i);
  const Point3D &operator[](size_t i) const;
  const PointVector &points() const;
  const std::list<Triangle> &triangles() const;
  static TriMesh readOBJ(std::string filename);
  void writeOBJ(std::string filename) const;

  const Triangle &closestTriangle(const Point3D &p) const;
  Point3D projectToTriangle(const Point3D &p, const Triangle &tri) const;

private:
  PointVector points_;
  std::list<Triangle> triangles_;
};

} // namespace Geometry
