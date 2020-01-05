#include <algorithm>

#include "geometry.hh"

namespace Geometry {

BSBasis::BSBasis() {
}

BSBasis::BSBasis(size_t degree, const DoubleVector &knots) : p_(degree), knots_(knots) {
}

size_t
BSBasis::degree() const {
  return p_;
}

void
BSBasis::setDegree(size_t degree) {
  p_ = degree;
}

const DoubleVector &
BSBasis::knots() const {
  return knots_;
}

DoubleVector &
BSBasis::knots() {
  return knots_;
}

double
BSBasis::low() const {
  return knots_[p_];
}

double
BSBasis::high() const {
  return knots_[knots_.size()-p_-1];
}

void
BSBasis::reverse() {
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
}

void
BSBasis::normalize() {
  size_t k = knots_.size();
  double low = knots_.front(), high = knots_.back(), len = high - low;
  for (size_t i = 0; i < k; ++i) {
    knots_[i] = (knots_[i] - low) / len;
  }
}

size_t
BSBasis::findSpan(double u) const {
  if(u >= knots_[knots_.size()-p_-1])
    return knots_.size() - p_ - 2;
  return (std::upper_bound(knots_.begin() + p_ + 1, knots_.end(), u) - knots_.begin()) - 1;
}

size_t
BSBasis::findSpanWithMultiplicity(double u, size_t &multi) const
{
  auto range = std::equal_range(knots_.begin(), knots_.end(), u);
  multi = range.second - range.first;

  if (u == knots_[knots_.size()-p_-1])
    return knots_.size() - p_ - 2;
  return (range.second - knots_.begin()) - 1;
}


void
BSBasis::basisFunctions(size_t i, double u, DoubleVector &coeff) const {
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
BSBasis::basisFunctionsAll(size_t i, double u, std::vector<DoubleVector> &coeff) const {
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

BSCurve::BSCurve() {
}

BSCurve::BSCurve(const PointVector &cpts) : n_(cpts.size() - 1), cp_(cpts) {
  size_t order = cpts.size();
  DoubleVector knots;
  basis_.setDegree(order - 1);
  basis_.knots().reserve(2 * order);
  basis_.knots().insert(basis_.knots().end(), order, 0.0);
  basis_.knots().insert(basis_.knots().end(), order, 1.0);
}

BSCurve::BSCurve(size_t degree, const DoubleVector &knots, const PointVector &cpts)
  : n_(cpts.size() - 1), cp_(cpts)
{
  basis_.setDegree(degree);
  basis_.knots() = knots;
}


void
BSCurve::derivativeControlPoints(size_t d, size_t r1, size_t r2,
                                 std::vector<PointVector> &dcp) const {
  size_t p = basis_.degree();
  const auto &knots = basis_.knots();
  dcp.clear(); dcp.resize(d + 1);
  size_t r = r2 - r1;
  dcp[0].reserve(r + 1);
  for(size_t i = 0; i <= r; ++i)
    dcp[0].push_back(cp_[r1+i]);
  for(size_t k = 1; k <= d; ++k) {
    dcp[k].reserve(r + 1 - k);
    size_t tmp = p - k + 1;
    for(size_t i = 0; i <= r - k; ++i)
      dcp[k].push_back((dcp[k-1][i+1] - dcp[k-1][i]) * tmp / (knots[r1+i+p+1] - knots[r1+i+k]));
  }
}


Point3D
BSCurve::eval(double u) const {
  size_t p = basis_.degree();
  size_t span = basis_.findSpan(u);
  DoubleVector coeff; basis_.basisFunctions(span, u, coeff);
  Point3D point(0.0, 0.0, 0.0);
  for(size_t i = 0; i <= p; ++i)
    point += cp_[span - p + i] * coeff[i];
  return point;
}

Point3D
BSCurve::eval(double u, size_t nr_der, VectorVector &der) const {
  size_t p = basis_.degree();
  size_t du = std::min(nr_der, p);
  der.clear();
  size_t span = basis_.findSpan(u);
  std::vector<DoubleVector> coeff; basis_.basisFunctionsAll(span, u, coeff);
  std::vector<PointVector> dcp; derivativeControlPoints(du, span - p, span, dcp);
  for(size_t k = 0; k <= du; ++k) {
    der.push_back(Vector3D(0.0, 0.0, 0.0));
    for(size_t j = 0; j <= p - k; ++j)
      der[k] += dcp[k][j] * coeff[p-k][j];
  }
  for(size_t k = p + 1; k <= nr_der; ++k)
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

const BSBasis &
BSCurve::basis() const {
  return basis_;
}

void
BSCurve::reverse() {
  basis_.reverse();
  std::reverse(cp_.begin(), cp_.end());
}

void
BSCurve::normalize() {
  basis_.normalize();
}

double
BSCurve::arcLength(double from, double to) const {
  size_t p = basis_.degree();
  const auto &knots = basis_.knots();
  if (to > knots[knots.size() - p - 1])
    to = knots[knots.size() - p - 1];
  if (from >= to)
    return 0.0;
  double next = std::min(to, knots[basis_.findSpan(from) + 1]);

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
  size_t p = basis_.degree();
  const auto &knots = basis_.knots();

  BSCurve result;
  result.basis_.setDegree(p); result.n_ = n_ + r;

  result.basis_.knots().reserve(knots.size() + r);
  std::copy_n(knots.begin(), k + 1, std::back_inserter(result.basis_.knots()));
  std::fill_n(std::back_inserter(result.basis_.knots()), r, u);
  std::copy(knots.begin() + k + 1, knots.end(), std::back_inserter(result.basis_.knots()));

  result.cp_.resize(cp_.size() + r);
  std::copy_n(cp_.begin(), k - p + 1, result.cp_.begin());
  std::copy(cp_.begin() + k - s, cp_.end(), result.cp_.begin() + r + k - s);

  PointVector tmp; tmp.reserve(p - s + 1);

  std::copy_n(cp_.begin() + k - p, p - s + 1, std::back_inserter(tmp));

  size_t L = k - p + 1;
  for (size_t j = 1; j <= r; ++j, ++L) {
    for (size_t i = 0; i <= p - j - s; ++i) {
      double alpha = (u - knots[L+i]) / (knots[i+k+1] - knots[L+i]);
      tmp[i] = tmp[i+1] * alpha + tmp[i] * (1.0 - alpha);
    }
    result.cp_[L] = tmp[0];
    result.cp_[k+r-j-s] = tmp[p-j-s];
  }
  std::copy_n(tmp.begin() + 1, p - s - 1 - r, result.cp_.begin() + L);

  return result;
}

BSCurve
BSCurve::insertKnot(double u, size_t r) const {
  size_t s;
  size_t k = basis_.findSpanWithMultiplicity(u, s);
  r = std::min(r, basis_.degree() - s);
  return insertKnot(u, k, s, r);
}

DoubleVector
BSCurve::intersectWithPlane(const Point3D &p, const Vector3D &n) const {
  DoubleVector result;
  BSCurve c = *this;
  for (size_t k = 1; k <= c.n_; ++k) {
    size_t deg = c.basis().degree();
    const auto &knots = c.basis().knots();
    double d1 = (c.cp_[k-1] - p) * n;
    double d2 = (c.cp_[k]   - p) * n;
    if (d1 * d2 <= 0.0) {
      if (d1 * d2 == 0.0 && d1 + d2 != 0.0 &&
          ((k == 1    && d1 == 0.0) ||
           (k == c.n_ && d2 == 0.0) ||
           (k < c.n_  && (c.cp_[k+1] - p) * n == 0)))
        continue;
      double k1 = knots[k], k2 = knots[k+deg];
      double g1 = k1, g2 = k2;
      for (size_t l = k + 1, le = k + deg; l != le; ++l) {
        g1 += knots[l];
        g2 += knots[l];
      }
      g1 /= deg; g2 /= deg;
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
        c = c.insertKnot(new_knot, deg - 1);
        if (c.n_ != old_size)
          --k;
      }
    }
  }
  return result;
}

} // namespace Geometry
