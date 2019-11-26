#include <algorithm>

#include "geometry.hh"

namespace Geometry {

BSCurve::BSCurve() {
}

BSCurve::BSCurve(const PointVector &cpts)
  : p_(cpts.size() - 1), n_(cpts.size() - 1), cp_(cpts) {
  size_t order = cpts.size();
  knots_.reserve(2 * order);
  knots_.insert(knots_.end(), order, 0.0);
  knots_.insert(knots_.end(), order, 1.0);
}

BSCurve::BSCurve(size_t degree, const DoubleVector &knots, const PointVector &cpts)
  : p_(degree), n_(cpts.size() - 1), knots_(knots), cp_(cpts) {
}

size_t
BSCurve::findSpan(double u) const {
  if(u == knots_[n_+1])
    return n_;
  return (std::upper_bound(knots_.begin() + p_ + 1, knots_.end(), u) - knots_.begin()) - 1;
}

size_t
BSCurve::findSpanWithMultiplicity(double u, size_t &multi) const
{
  auto range = std::equal_range(knots_.begin(), knots_.end(), u);
  multi = range.second - range.first;

  if (u == knots_[n_+1])
    return n_;
  return (range.second - knots_.begin()) - 1;
}


void
BSCurve::basisFunctions(size_t i, double u, DoubleVector &coeff) const {
  coeff.clear(); coeff.reserve(p_ + 1);
  coeff.push_back(1.0);
  DoubleVector left(p_ + 1), right(p_ + 1);
  for(size_t j = 1; j <= p_; ++j) {
    left[j]  = u - knots_[i+1-j];
    right[j] = knots_[i+j] - u;
    double saved = 0.0;
    for(size_t r = 0; r < j; ++r) {
      double tmp = coeff[r] / (right[r+1] + left[j-r]);
      coeff[r] = saved + tmp * right[r+1];
      saved = tmp * left[j-r];
    }
    coeff.push_back(saved);
  }
}

void
BSCurve::derivativeControlPoints(size_t d, size_t r1, size_t r2,
                                 std::vector<PointVector> &dcp) const {
  dcp.clear(); dcp.resize(d + 1);
  size_t r = r2 - r1;
  dcp[0].reserve(r + 1);
  for(size_t i = 0; i <= r; ++i)
    dcp[0].push_back(cp_[r1+i]);
  for(size_t k = 1; k <= d; ++k) {
    dcp[k].reserve(r + 1 - k);
    size_t tmp = p_ - k + 1;
    for(size_t i = 0; i <= r - k; ++i)
      dcp[k].push_back((dcp[k-1][i+1] - dcp[k-1][i]) * tmp / (knots_[r1+i+p_+1] - knots_[r1+i+k]));
  }
}

void
BSCurve::basisFunctionsAll(size_t i, double u, std::vector<DoubleVector> &coeff) const {
  coeff.clear(); coeff.resize(p_ + 1);
  coeff[0].push_back(1.0);
  DoubleVector left(p_ + 1), right(p_ + 1);
  for(size_t j = 1; j <= p_; ++j) {
    coeff[j].reserve(j + 1);
    left[j]  = u - knots_[i+1-j];
    right[j] = knots_[i+j] - u;
    double saved = 0.0;
    for(size_t r = 0; r < j; ++r) {
      double tmp = coeff[j-1][r] / (right[r+1] + left[j-r]);
      coeff[j].push_back(saved + tmp * right[r+1]);
      saved = tmp * left[j-r];
    }
    coeff[j].push_back(saved);
  }
}

Point3D
BSCurve::eval(double u) const {
  size_t span = findSpan(u);
  DoubleVector coeff; basisFunctions(span, u, coeff);
  Point3D point(0.0, 0.0, 0.0);
  for(size_t i = 0; i <= p_; ++i)
    point += cp_[span - p_ + i] * coeff[i];
  return point;
}

Point3D
BSCurve::eval(double u, size_t nr_der, VectorVector &der) const {
  size_t du = std::min(nr_der, p_);
  der.clear();
  size_t span = findSpan(u);
  std::vector<DoubleVector> coeff; basisFunctionsAll(span, u, coeff);
  std::vector<PointVector> dcp; derivativeControlPoints(du, span - p_, span, dcp);
  for(size_t k = 0; k <= du; ++k) {
    der.push_back(Vector3D(0.0, 0.0, 0.0));
    for(size_t j = 0; j <= p_ - k; ++j)
      der[k] += dcp[k][j] * coeff[p_-k][j];
  }
  for(size_t k = p_ + 1; k <= nr_der; ++k)
    der.push_back(Vector3D(0.0, 0.0, 0.0));
  return der[0];
}

const PointVector &
BSCurve::controlPoints() const {
  return cp_;
}

PointVector &
BSCurve::controlPoints() {
  return cp_;
}

void
BSCurve::reverse() {
  size_t k = knots_.size();
  DoubleVector new_knots;
  new_knots.reserve(k);

  double curr = knots_.front();
  for (size_t i = 1, j = k - 1; i < k; ++i, --j) {
    new_knots.push_back(curr);
    curr += knots_[j] - knots_[j-1];
  }
  new_knots.push_back(curr);

  knots_ = new_knots;
  std::reverse(cp_.begin(), cp_.end());
}

void
BSCurve::normalize() {
  size_t k = knots_.size();
  double low = knots_.front(), high = knots_.back(), len = high - low;
  for (size_t i = 0; i < k; ++i) {
    knots_[i] = (knots_[i] - low) / len;
  }
}

double
BSCurve::arcLength(double from, double to) const {
  if (to > knots_[knots_.size() - p_ - 1])
    to = knots_[knots_.size() - p_ - 1];
  if (from >= to)
    return 0.0;
  double next = std::min(to, knots_[findSpan(from) + 1]);

  // Estimate each knot interval using Gaussian quadratures
  const static double gauss[] = {-0.861136312, 0.347854845,
                                 -0.339981044, 0.652145155,
                                 0.339981044, 0.652145155,
                                 0.861136312, 0.347854845};

  double sum = 0.0;
  for (size_t i = 0; i < 8; i += 2) {
    VectorVector der;
    double u = ((next - from) * gauss[i] + from + next) * 0.5;
    eval(u, 1, der);
    sum += der[1].norm() * gauss[i+1] * (next - from) * 0.5;
  }

  return sum + arcLength(next, to);
}

BSCurve
BSCurve::insertKnot(double u, size_t k, size_t s, size_t r) const {
  BSCurve result;
  result.p_ = p_; result.n_ = n_ + r;

  result.knots_.reserve(knots_.size() + r);
  std::copy_n(knots_.begin(), k + 1, std::back_inserter(result.knots_));
  std::fill_n(std::back_inserter(result.knots_), r, u);
  std::copy(knots_.begin() + k + 1, knots_.end(), std::back_inserter(result.knots_));

  result.cp_.resize(cp_.size() + r);
  std::copy_n(cp_.begin(), k - p_ + 1, result.cp_.begin());
  std::copy(cp_.begin() + k - s, cp_.end(), result.cp_.begin() + r + k - s);

  PointVector tmp; tmp.reserve(p_ - s + 1);

  std::copy_n(cp_.begin() + k - p_, p_ - s + 1, std::back_inserter(tmp));

  size_t L = k - p_ + 1;
  for (size_t j = 1; j <= r; ++j, ++L) {
    for (size_t i = 0; i <= p_ - j - s; ++i) {
      double alpha = (u - knots_[L+i]) / (knots_[i+k+1] - knots_[L+i]);
      tmp[i] = tmp[i+1] * alpha + tmp[i] * (1.0 - alpha);
    }
    result.cp_[L] = tmp[0];
    result.cp_[k+r-j-s] = tmp[p_-j-s];
  }
  std::copy_n(tmp.begin() + 1, p_ - s - 1 - r, result.cp_.begin() + L);

  return result;
}

BSCurve
BSCurve::insertKnot(double u, size_t r) const {
  size_t s;
  size_t k = findSpanWithMultiplicity(u, s);
  r = std::min(r, p_ - s);
  return insertKnot(u, k, s, r);
}

DoubleVector
BSCurve::intersectWithPlane(const Point3D &p, const Vector3D &n) const {
  DoubleVector result;
  BSCurve c = *this;
  for (size_t k = 1; k <= c.n_; ++k) {
    double d1 = (c.cp_[k-1] - p) * n;
    double d2 = (c.cp_[k]   - p) * n;
    if (d1 * d2 <= 0.0) {
      if (d1 * d2 == 0.0 && d1 + d2 != 0.0 &&
          ((k == 1    && d1 == 0.0) ||
           (k == c.n_ && d2 == 0.0) ||
           (k < c.n_  && (c.cp_[k+1] - p) * n == 0)))
        continue;
      double k1 = c.knots_[k], k2 = c.knots_[k+c.p_];
      double g1 = k1, g2 = k2;
      for (size_t l = k + 1, le = k + c.p_; l != le; ++l) {
        g1 += c.knots_[l];
        g2 += c.knots_[l];
      }
      g1 /= c.p_; g2 /= c.p_;
      double ratio = d1 == 0.0 && d2 == 0.0
        ? (g1 == k1 ? 0.0 : (g2 == k2 ? 1.0 : 0.5))
        : std::abs(d1 / (d1 - d2));
      double new_knot = g1 * (1.0 - ratio) + g2 * ratio;
      if (std::abs((c.eval(new_knot) - p) * n) < epsilon) {
        result.push_back(new_knot);
        for (; k <= c.n_; ++k) {
          double dd = (c.cp_[k] - p) * n;
          if (dd != 0.0)
            break;
        }
      } else {
        size_t old_size = c.n_;
        c = c.insertKnot(new_knot, c.p_ - 1);
        if (c.n_ != old_size)
          --k;
      }
    }
  }
  return result;
}

} // namespace Geometry
