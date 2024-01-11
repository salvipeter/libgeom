#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "geometry.hh"

namespace Geometry {

Matrix2x2::Matrix2x2(const double *values) {
  std::copy_n(values, 4, m_.begin());
}

Matrix2x2::Matrix2x2(std::initializer_list<double> values) {
  if (values.size() != 4)
    throw std::invalid_argument("Matrix2x2 needs 4 values");
  std::copy(values.begin(), values.end(), m_.begin());
}

Matrix2x2
Matrix2x2::identity() {
  Matrix2x2 I;
  I.m_[0] = 1.0;
  I.m_[1] = 0.0;
  I.m_[2] = 0.0;
  I.m_[3] = 1.0;
  return I;
}

Matrix2x2
Matrix2x2::rotation(double angle) {
  double c = std::cos(angle), s = std::sin(angle);
  Matrix2x2 A;
  A.m_[0] = c;
  A.m_[1] = -s;
  A.m_[2] = s;
  A.m_[3] = c;
  return A;
}

Matrix2x2
Matrix2x2::operator+(const Matrix2x2 &m) const {
  Matrix2x2 result;
  for (size_t i = 0; i < 4; ++i)
    result.m_[i] = m_[i] + m.m_[i];
  return result;
}

Matrix2x2 &
Matrix2x2::operator+=(const Matrix2x2 &m) {
  for (size_t i = 0; i < 4; ++i)
    m_[i] += m.m_[i];
  return *this;
}

Matrix2x2
Matrix2x2::operator*(double x) const {
  Matrix2x2 result;
  for (size_t i = 0; i < 4; ++i)
    result.m_[i] = m_[i] * x;
  return result;
}

Matrix2x2 &
Matrix2x2::operator*=(double x) {
  for (size_t i = 0; i < 4; ++i)
    m_[i] *= x;
  return *this;
}

Vector2D
Matrix2x2::operator*(const Vector2D &v) const {
  return {
    m_[0] * v[0] + m_[1] * v[1],
    m_[2] * v[0] + m_[3] * v[1]
  };
}

Matrix2x2
Matrix2x2::operator*(const Matrix2x2 &m) const {
  Matrix2x2 result;
  result.m_[0] = m_[0] * m.m_[0] + m_[1] * m.m_[2];
  result.m_[1] = m_[0] * m.m_[1] + m_[1] * m.m_[3];
  result.m_[2] = m_[2] * m.m_[0] + m_[3] * m.m_[2];
  result.m_[3] = m_[2] * m.m_[1] + m_[3] * m.m_[3];
  return result;
}

Matrix2x2 &
Matrix2x2::operator*=(const Matrix2x2 &m) {
  *this = (*this) * m;
  return *this;
}

double Matrix2x2::trace() const {
  return m_[0] + m_[3];
}

double Matrix2x2::determinant() const {
  return m_[0] * m_[3] - m_[1] * m_[2];
}

Matrix2x2 Matrix2x2::adjugate() const {
  Matrix2x2 r;
  r.m_[0] = m_[3];
  r.m_[1] = -m_[1];
  r.m_[2] = -m_[2];
  r.m_[3] = m_[0];
  return r;
}

Matrix2x2 Matrix2x2::transpose() const {
  Matrix2x2 r;
  r.m_[0] = m_[0];
  r.m_[1] = m_[2];
  r.m_[2] = m_[1];
  r.m_[3] = m_[3];
  return r;
}

Matrix2x2 Matrix2x2::inverse() const {
  return adjugate() * (1 / determinant());
}

const double &Matrix2x2::operator()(size_t i, size_t j) const {
  return m_[i * 2 + j];
}

double &Matrix2x2::operator()(size_t i, size_t j) {
  return m_[i * 2 + j];
}

const double *Matrix2x2::data() const {
  return &m_[0];
}


} // namespace Geometry
