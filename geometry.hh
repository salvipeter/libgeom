#pragma once

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
using DoubleMatrix = std::vector<DoubleVector>;
using Vector2DVector = std::vector<Vector2D>;
using Vector2DMatrix = std::vector<Vector2DVector>;
using VectorVector = std::vector<Vector3D>;
using VectorMatrix = std::vector<VectorVector>;
using Point2DVector = std::vector<Point2D>;
using Point2DMatrix = std::vector<Point2DVector>;
using PointVector = std::vector<Point3D>;
using PointMatrix = std::vector<PointVector>;

class Vector2D {
public:
  // Constructors
  Vector2D();
  explicit Vector2D(const double *coords);
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
  Vector2D normalized() const;

private:
  std::array<double, 2> v_;
};

std::ostream &operator<<(std::ostream &os, const Vector2D &v);
std::istream &operator>>(std::istream &is, Vector2D &v);

class Vector3D {
public:
  // Constructors
  Vector3D();
  explicit Vector3D(const double *coords);
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
  Vector3D normalized() const;

private:
  std::array<double, 3> v_;
};

std::ostream &operator<<(std::ostream &os, const Vector3D &v);
std::istream &operator>>(std::istream &is, Vector3D &v);

class Matrix3x3 {
public:
  // Constructors
  Matrix3x3() = default;
  explicit Matrix3x3(const double *values);        // row-major
  Matrix3x3(std::initializer_list<double> values); // row-major

  // Special matrices
  static Matrix3x3 identity();
  static Matrix3x3 rotation(const Vector3D &axis, double angle);

  const double &operator()(size_t i, size_t j) const;
  double &operator()(size_t i, size_t j);
  const double *data() const;

  // Arithmetic
  Matrix3x3 operator+(const Matrix3x3 &m) const;
  Matrix3x3 &operator+=(const Matrix3x3 &m);
  Matrix3x3 operator*(double x) const;
  Matrix3x3 &operator*=(double x);
  Vector3D operator*(const Vector3D &v) const;
  Matrix3x3 operator*(const Matrix3x3 &m) const;
  Matrix3x3 &operator*=(const Matrix3x3 &m);

  double trace() const;
  Matrix3x3 adjugate() const;
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
  double low() const;
  double high() const;

  // Utilities
  void reverse();
  void normalize();
  size_t findSpan(double u) const;
  size_t findSpanWithMultiplicity(double u, size_t &multi) const;
  void basisFunctions(size_t i, double u, DoubleVector &coeff) const;
  void basisFunctionsAll(size_t i, double u, DoubleMatrix &coeff) const;
  void basisFunctionDerivatives(size_t i, double u, size_t nr_der, DoubleMatrix &coeff) const;

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
  BSCurve insertKnot(double u, size_t k, size_t s, size_t r) const;

  size_t n_;
  BSBasis basis_;
  PointVector cp_;
};

class BSSurface {
public:
  // Constructors
  BSSurface();
  BSSurface(const BSSurface &) = default;
  BSSurface(size_t deg_u, size_t deg_v, const PointVector &cpts);
  BSSurface(size_t deg_u, size_t deg_v, const DoubleVector &knots_u, const DoubleVector &knots_v,
            const PointVector &cpts);
  BSSurface &operator=(const BSSurface &) = default;

  // Evaluation
  Point3D eval(double u, double v) const;
  Point3D eval(double u, double v, size_t nr_der, VectorMatrix &der) const;

  // Coordinates
  std::array<size_t, 2> numControlPoints() const;
  Point3D controlPoint(size_t i, size_t j) const;
  Point3D &controlPoint(size_t i, size_t j);
  const PointVector &controlPoints() const;
  PointVector &controlPoints();

  // Parameterization
  const BSBasis &basisU() const;
  const BSBasis &basisV() const;
  void swapUV();
  void reverseU();
  void reverseV();
  void normalize();

  // Algorithms
  BSSurface insertKnotU(double u, size_t r) const;
  BSSurface insertKnotV(double v, size_t r) const;

private:
  BSSurface insertKnotU(double u, size_t k, size_t s, size_t r) const;
  BSSurface insertKnotV(double v, size_t k, size_t s, size_t r) const;

  size_t n_u_, n_v_;
  BSBasis basis_u_, basis_v_;
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
  TriMesh &insert(const TriMesh &other, double tolerance);

  // I/O
  Point3D &operator[](size_t i);
  const Point3D &operator[](size_t i) const;
  const PointVector &points() const;
  const std::list<Triangle> &triangles() const;
  static TriMesh readOBJ(std::string filename);
  void writeOBJ(std::string filename) const;
  void writeSTL(std::string filename) const;

  const Triangle &closestTriangle(const Point3D &p) const;
  Point3D projectToTriangle(const Point3D &p, const Triangle &tri) const;

private:
  PointVector points_;
  std::list<Triangle> triangles_;
};

} // namespace Geometry
