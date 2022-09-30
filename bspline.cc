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

  if (u >= knots_[knots_.size()-p_-1])
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
BSBasis::basisFunctionsAll(size_t i, double u, DoubleMatrix &coeff) const {
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

void
BSBasis::basisFunctionDerivatives(size_t i, double u, size_t nr_der, DoubleMatrix &coeff) const {
  coeff.clear(); coeff.resize(nr_der + 1);
  DoubleMatrix ndu(p_ + 1);
  ndu[0].push_back(1);
  DoubleVector left(p_ + 1), right(p_ + 1);
  for (size_t j = 1; j <= p_; ++j) {
    left[j] = u - knots_[i+1-j];
    right[j] = knots_[i+j] - u;
    double saved = 0;
    for (size_t r = 0; r < j; ++r) {
      ndu[j].push_back(right[r+1] + left[j-r]);
      double tmp = ndu[r][j-1] / ndu[j][r];
      ndu[r].push_back(saved + right[r+1] * tmp);
      saved = tmp * left[j-r];
    }
    ndu[j].push_back(saved);
  }
  for (size_t j = 0; j <= p_; ++j)
    coeff[0].push_back(ndu[j][p_]);
  DoubleVector a[2]; a[0].resize(p_ + 1); a[1].resize(p_ + 1);
  for (size_t r = 0; r <= p_; ++r) {
    size_t s1 = 0, s2 = 1;
    a[0][0] = 1;
    for (size_t k = 1; k <= nr_der; ++k) {
      double d = 0;
      size_t pk = p_ - k;
      if (r >= k) {
        a[s2][0] = a[s1][0] / ndu[pk+1][r-k];
        d = a[s2][0] * ndu[r-k][pk];
      }
      size_t j1 = r >= k - 1 ? 1 : k - r;
      size_t j2 = r <= pk + 1 ? k - 1 : p_ - r;
      for (size_t j = j1; j <= j2; ++j) {
        a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu[pk+1][r+j-k];
        d += a[s2][j] * ndu[r+j-k][pk];
      }
      if (r <= pk) {
        a[s2][k] = -a[s1][k-1] / ndu[pk+1][r];
        d += a[s2][k] * ndu[r][pk];
      }
      coeff[k].push_back(d);
      std::swap(s1, s2);
    }
  }
  size_t r = p_;
  for (size_t k = 1; k <= nr_der; ++k) {
    for (size_t j = 0; j <= p_; ++j)
      coeff[k][j] *= r;
    r *= p_ - k;
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
  std::vector<DoubleVector> coeff; basis_.basisFunctionDerivatives(span, u, du, coeff);
  for(size_t k = 0; k <= du; ++k) {
    der.emplace_back(0.0, 0.0, 0.0);
    for(size_t j = 0; j <= p; ++j)
      der[k] += cp_[span-p+j] * coeff[k][j];
  }
  for(size_t k = p + 1; k <= nr_der; ++k)
    der.emplace_back(0.0, 0.0, 0.0);
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
  if (p > s + r + 1)
    std::copy_n(tmp.begin() + 1, p - s - 1 - r, result.cp_.begin() + L);

  return result;
}

BSCurve
BSCurve::insertKnot(double u, size_t r) const {
  size_t s;
  size_t k = basis_.findSpanWithMultiplicity(u, s);
  if (s >= basis_.degree())
    return *this;
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

BSSurface::BSSurface() {
}

BSSurface::BSSurface(size_t deg_u, size_t deg_v, const PointVector &cpts)
  : n_u_(deg_u), n_v_(deg_v), cp_(cpts)
{
  basis_u_.setDegree(deg_u);
  basis_v_.setDegree(deg_v);
  basis_u_.knots().reserve(2 * (deg_u + 1));
  basis_u_.knots().insert(basis_u_.knots().end(), deg_u + 1, 0.0);
  basis_u_.knots().insert(basis_u_.knots().end(), deg_u + 1, 1.0);
  basis_v_.knots().reserve(2 * (deg_v + 1));
  basis_v_.knots().insert(basis_v_.knots().end(), deg_v + 1, 0.0);
  basis_v_.knots().insert(basis_v_.knots().end(), deg_v + 1, 1.0);
}

BSSurface::BSSurface(size_t deg_u, size_t deg_v,
                     const DoubleVector &knots_u, const DoubleVector &knots_v,
                     const PointVector &cpts)
  : n_u_(knots_u.size() - deg_u - 2), n_v_(knots_v.size() - deg_v - 2), cp_(cpts)
{
  basis_u_.setDegree(deg_u);
  basis_v_.setDegree(deg_v);
  basis_u_.knots() = knots_u;
  basis_v_.knots() = knots_v;
}

Point3D
BSSurface::eval(double u, double v) const {
  size_t p_u = basis_u_.degree(), p_v = basis_v_.degree();
  size_t span_u = basis_u_.findSpan(u), span_v = basis_v_.findSpan(v);
  DoubleVector coeff_u, coeff_v;
  basis_u_.basisFunctions(span_u, u, coeff_u);
  basis_v_.basisFunctions(span_v, v, coeff_v);
  Point3D point(0.0, 0.0, 0.0);
  for (size_t i = 0; i <= p_u; ++i) {
    size_t base = (span_u - p_u + i) * (n_v_ + 1);
    for (size_t j = 0; j <= p_v; ++j)
      point += cp_[base+span_v-p_v+j] * coeff_u[i] * coeff_v[j];
  }
  return point;
}

Point3D
BSSurface::eval(double u, double v, size_t nr_der, VectorMatrix &der) const {
  der.clear(); der.resize(nr_der + 1);
  size_t p_u = basis_u_.degree(), p_v = basis_v_.degree();
  size_t du = std::min(nr_der, p_u), dv = std::min(nr_der, p_v);
  size_t span_u = basis_u_.findSpan(u), span_v = basis_v_.findSpan(v);
  DoubleMatrix coeff_u, coeff_v;
  basis_u_.basisFunctionDerivatives(span_u, u, du, coeff_u);
  basis_v_.basisFunctionDerivatives(span_v, v, dv, coeff_v);
  VectorVector tmp(p_v + 1);
  for (size_t k = 0; k <= du; ++k) {
    for (size_t s = 0; s <= p_v; ++s) {
      tmp[s] = Vector3D(0.0, 0.0, 0.0);
      for (size_t r = 0; r <= p_u; ++r)
        tmp[s] += cp_[(span_u-p_u+r)*(n_v_+1)+span_v-p_v+s] * coeff_u[k][r];
    }
    size_t dd = std::min(nr_der - k, dv);
    for (size_t l = 0; l <= dd; ++l) {
      Vector3D point(0.0, 0.0, 0.0);
      for (size_t s = 0; s <= p_v; ++s)
        point += tmp[s] * coeff_v[l][s];
      der[k].push_back(point);
    }
  }
  for (size_t k = p_u + 1; k <= nr_der; ++k)
    for (size_t l = 0; l <= nr_der - k; ++l)
      der[k].emplace_back(0.0, 0.0, 0.0);
  for (size_t l = p_v + 1; l <= nr_der; ++l)
    for (size_t k = 0; k <= nr_der - l; ++k)
      der[k].emplace_back(0.0, 0.0, 0.0);
  return der[0][0];
}

const PointVector &
BSSurface::controlPoints() const {
  return cp_;
}

std::array<size_t, 2>
BSSurface::numControlPoints() const {
  return { n_u_ + 1, n_v_ + 1 };
}

Point3D
BSSurface::controlPoint(size_t i, size_t j) const {
  return cp_[i*(n_v_+1)+j];
}

Point3D &
BSSurface::controlPoint(size_t i, size_t j) {
  return cp_[i*(n_v_+1)+j];
}

PointVector &
BSSurface::controlPoints() {
  return cp_;
}

const BSBasis &
BSSurface::basisU() const {
  return basis_u_;
}

const BSBasis &
BSSurface::basisV() const {
  return basis_v_;
}

void
BSSurface::swapUV() {
  std::swap(n_u_, n_v_);
  std::swap(basis_u_, basis_v_);
  auto cp_old = cp_;
  for (size_t i = 0, index = 0; i <= n_u_; ++i)
    for (size_t j = 0; j <= n_v_; ++j, ++index)
      cp_[index] = cp_old[j*(n_u_+1)+i];
}

void
BSSurface::reverseU() {
  basis_u_.reverse();
  for (size_t i = 0; i < (n_u_ + 1) / 2; ++i) {
    size_t i2 = n_u_ - i;
    for (size_t j = 0; j <= n_v_; ++j)
      std::swap(cp_[i*(n_v_+1)+j], cp_[i2*(n_v_+1)+j]);
  }
}

void BSSurface::reverseV() {
  basis_v_.reverse();
  for (size_t i = 0; i <= n_u_; ++i) {
    size_t base = i * (n_v_ + 1);
    for (size_t j = 0; j < (n_v_ + 1) / 2; ++j)
      std::swap(cp_[base+j], cp_[base+n_v_-j]);
  }
}

void BSSurface::normalize() {
  basis_u_.normalize();
  basis_v_.normalize();
}

BSSurface
BSSurface::insertKnotU(double u, size_t k, size_t s, size_t r) const {
  size_t p = basis_u_.degree();
  const auto &knots = basis_u_.knots();

  BSSurface result;
  result.basis_v_ = basis_v_;
  result.n_v_ = n_v_;

  result.basis_u_.setDegree(p); result.n_u_ = n_u_ + r;

  result.basis_u_.knots().reserve(knots.size() + r);
  std::copy_n(knots.begin(), k + 1, std::back_inserter(result.basis_u_.knots()));
  std::fill_n(std::back_inserter(result.basis_u_.knots()), r, u);
  std::copy(knots.begin() + k + 1, knots.end(), std::back_inserter(result.basis_u_.knots()));

  size_t ncol = n_v_ + 1;
  result.cp_.resize((n_u_ + r + 1) * ncol);
  std::copy_n(cp_.begin(), (k - p + 1) * ncol, result.cp_.begin());
  std::copy(cp_.begin() + (k - s) * ncol, cp_.end(), result.cp_.begin() + (r + k - s) * ncol);

  PointVector tmp; tmp.reserve((p - s + 1) * ncol);

  std::copy_n(cp_.begin() + (k - p) * ncol, (p - s + 1) * ncol, std::back_inserter(tmp));

  size_t L = k - p + 1;
  for (size_t j = 1; j <= r; ++j, ++L) {
    for (size_t i = 0; i <= p - j - s; ++i) {
      double alpha = (u - knots[L+i]) / (knots[i+k+1] - knots[L+i]);
      for (size_t l = 0; l < ncol; ++l)
        tmp[i*ncol+l] = tmp[(i+1)*ncol+l] * alpha + tmp[i*ncol+l] * (1.0 - alpha);
    }
    for (size_t l = 0; l < ncol; ++l) {
      result.cp_[L*ncol+l] = tmp[l];
      result.cp_[(k+r-j-s)*ncol+l] = tmp[(p-j-s)*ncol+l];
    }
  }
  if (p > s + r + 1)
    std::copy_n(tmp.begin() + ncol, (p - s - 1 - r) * ncol, result.cp_.begin() + L * ncol);

  return result;
}

BSSurface
BSSurface::insertKnotV(double v, size_t k, size_t s, size_t r) const {
  auto copy = *this;
  copy.swapUV();
  auto result = copy.insertKnotU(v, k, s, r);
  result.swapUV();
  return result;
}

BSSurface
BSSurface::insertKnotU(double u, size_t r) const {
  size_t s;
  size_t k = basis_u_.findSpanWithMultiplicity(u, s);
  if (s >= basis_u_.degree())
    return *this;
  r = std::min(r, basis_u_.degree() - s);
  return insertKnotU(u, k, s, r);
}

BSSurface
BSSurface::insertKnotV(double v, size_t r) const {
  size_t s;
  size_t k = basis_v_.findSpanWithMultiplicity(v, s);
  if (s >= basis_v_.degree())
    return *this;
  r = std::min(r, basis_v_.degree() - s);
  return insertKnotV(v, k, s, r);
}

} // namespace Geometry
