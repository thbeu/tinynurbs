#include <tinynurbs/tinynurbs.h>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <cmath>

#include "catch.hpp"

using namespace std;

tinynurbs::RationalSurface3f getHemisphere() {
    tinynurbs::RationalSurface3f srf;
    srf.degree_u = 3;
    srf.degree_v = 3;
    srf.knots_u = {0, 0, 0, 0, 1, 1, 1, 1};
    srf.knots_v = {0, 0, 0, 0, 1, 1, 1, 1};
    // 4x4 grid (tinynurbs::array2) of control points and weights
    // https://www.geometrictools.com/Documentation/NURBSCircleSphere.pdf
    srf.control_points = {4, 4, 
                          {glm::vec3(0, 0, 1), glm::vec3(0, 0, 1), glm::vec3(0, 0, 1), glm::vec3(0, 0, 1),
                           glm::vec3(2, 0, 1), glm::vec3(2, 4, 1),  glm::vec3(-2, 4, 1),  glm::vec3(-2, 0, 1),
                           glm::vec3(2, 0, -1), glm::vec3(2, 4, -1), glm::vec3(-2, 4, -1), glm::vec3(-2, 0, -1),
                           glm::vec3(0, 0, -1), glm::vec3(0, 0, -1), glm::vec3(0, 0, -1), glm::vec3(0, 0, -1)
                          }
                         };
    srf.weights = {4, 4,
                   {1,       1.f/3.f, 1.f/3.f, 1,
                    1.f/3.f, 1.f/9.f, 1.f/9.f, 1.f/3.f,
                    1.f/3.f, 1.f/9.f, 1.f/9.f, 1.f/3.f,
                    1,       1.f/3.f, 1.f/3.f, 1
                   }
                 };
    return srf;
}

TEST_CASE("surfacePoint (rational)", "[surface, rational, evaluate]")
{
    auto srf = getHemisphere();
    glm::vec3 pt1 = tinynurbs::surfacePoint(srf, 0.f, 0.f);
    glm::vec3 pt2 = tinynurbs::surfacePoint(srf, 1.f, 1.f);
    REQUIRE(pt1.x == Approx(0));
    REQUIRE(pt1.y == Approx(0));
    REQUIRE(pt1.z == Approx(1));
    REQUIRE(pt2.x == Approx(0));
    REQUIRE(pt2.y == Approx(0));
    REQUIRE(pt2.z == Approx(-1));
}

TEST_CASE("surfaceTangent (rational)", "[surface, rational, evaluate]")
{
    auto srf = getHemisphere();
    glm::vec3 tgt_u, tgt_v;
    std::tie(tgt_u, tgt_v) = tinynurbs::surfaceTangent(srf, 0.f, 0.f);
    REQUIRE(glm::length(tgt_u) == Approx(1));
    REQUIRE(tgt_u.x == Approx(1));
    REQUIRE(tgt_u.y == Approx(0));
    REQUIRE(tgt_u.z == Approx(0));
    // tgt_v should be a zero vector since the pole of the hemisphere is pinched
    REQUIRE(glm::length(tgt_v) == Approx(0));
    REQUIRE(tgt_v.x == Approx(0));
    REQUIRE(tgt_v.y == Approx(0));
    REQUIRE(tgt_v.z == Approx(0));
}


TEST_CASE("surfaceNormal (rational)", "[surface, rational, evaluate]")
{
    auto srf = getHemisphere();
    glm::vec3 n = tinynurbs::surfaceNormal(srf, 0.5f, 0.5f);
    REQUIRE(glm::length(n) == Approx(1));
    REQUIRE(n.x == Approx(0));
    REQUIRE(n.y == Approx(-1));
    REQUIRE(n.z == Approx(0));

    // Normal is zero vector since the pole is pinched
    n = tinynurbs::surfaceNormal(srf, 0.f, 0.f);
    REQUIRE(glm::length(n) == Approx(0));
    REQUIRE(n.x == Approx(0));
    REQUIRE(n.y == Approx(0));
    REQUIRE(n.z == Approx(0));
}

TEST_CASE("surfaceIsValid (rational)", "[surface, rational, check]")
{
    auto srf = getHemisphere();
    bool is_valid = tinynurbs::surfaceIsValid(srf);
    REQUIRE(is_valid == true);

    srf = getHemisphere();
    srf.degree_u = 5;
    is_valid = tinynurbs::surfaceIsValid(srf);
    REQUIRE(is_valid == false);

    srf = getHemisphere();
    srf.degree_v = 4;
    is_valid = tinynurbs::surfaceIsValid(srf);
    REQUIRE(is_valid == false);
}

TEST_CASE("surfaceInsertKnotU (rational)", "[surface, rational, modify]")
{
    auto srf = getHemisphere();
    unsigned int repeat = 2;

    glm::vec3 pt = tinynurbs::surfacePoint(srf, 0.25f, 0.5f);
    
    size_t n_knots_prev = srf.knots_u.size();
    size_t n_control_points_prev = srf.control_points.rows();

    auto new_srf = tinynurbs::surfaceKnotInsertU(srf, 0.25f, repeat);
    glm::vec3 new_pt = tinynurbs::surfacePoint(new_srf, 0.25f, 0.5f);

    size_t n_knots_curr = new_srf.knots_u.size();
    size_t n_control_points_curr = new_srf.control_points.rows();

    REQUIRE((n_knots_prev + repeat) == n_knots_curr);
    REQUIRE((n_control_points_prev + repeat) == n_control_points_curr);
    REQUIRE(pt.x == Approx(new_pt.x));
    REQUIRE(pt.y == Approx(new_pt.y));
    REQUIRE(pt.z == Approx(new_pt.z));
}

TEST_CASE("surfaceInsertKnotV (rational)", "[surface, rational, modify]")
{
    auto srf = getHemisphere();
    unsigned int repeat = 2;

    glm::vec3 pt = tinynurbs::surfacePoint(srf, 0.25f, 0.5f);
    
    size_t n_knots_prev = srf.knots_v.size();
    size_t n_control_points_prev = srf.control_points.cols();

    auto new_srf = tinynurbs::surfaceKnotInsertV(srf, 0.5f, repeat);
    glm::vec3 new_pt = tinynurbs::surfacePoint(new_srf, 0.25f, 0.5f);

    size_t n_knots_curr = new_srf.knots_v.size();
    size_t n_control_points_curr = new_srf.control_points.cols();

    REQUIRE((n_knots_prev + repeat) == n_knots_curr);
    REQUIRE((n_control_points_prev + repeat) == n_control_points_curr);
    REQUIRE(pt.x == Approx(new_pt.x));
    REQUIRE(pt.y == Approx(new_pt.y));
    REQUIRE(pt.z == Approx(new_pt.z));
}

TEST_CASE("surfaceSplitU (rational)", "[surface, rational, modify]")
{
    auto srf = getHemisphere();
    float u = 0.25f;
    tinynurbs::RationalSurface3f left, right;
    std::tie(left, right) = tinynurbs::surfaceSplitU(srf, u);

    bool is_valid_l = tinynurbs::surfaceIsValid(left);
    bool is_valid_r = tinynurbs::surfaceIsValid(right);

    REQUIRE(is_valid_l == true);
    REQUIRE(is_valid_r == true);

    REQUIRE(left.degree_u == srf.degree_u);
    REQUIRE(right.degree_u == srf.degree_u);
    REQUIRE(left.degree_v == srf.degree_v);
    REQUIRE(right.degree_v == srf.degree_v);

    for (int i = 0; i < left.degree_u + 1; ++i) {
        int d = left.knots_u.size() - (left.degree_u + 1);
        REQUIRE(left.knots_u[d+i] == Approx(u));
    }

    for (int i = 0; i < right.degree_u + 1; ++i) {
        REQUIRE(right.knots_u[i] == Approx(u));
    }

    tinynurbs::surfaceSaveOBJ("left_rational_u.obj", left);
    tinynurbs::surfaceSaveOBJ("right_rational_u.obj", right);

    glm::vec3 pt1 = tinynurbs::surfacePoint(srf, left.knots_u[left.knots_u.size() - 1], 0.f);
    glm::vec3 pt2 = tinynurbs::surfacePoint(srf, right.knots_u[0], 0.f);
    REQUIRE(pt1.x == Approx(pt2.x));
    REQUIRE(pt1.y == Approx(pt2.y));
    REQUIRE(pt1.z == Approx(pt2.z));
}

TEST_CASE("surfaceSplitV (rational)", "[surface, rational, modify]")
{
    auto srf = getHemisphere();
    float v = 0.25f;
    tinynurbs::RationalSurface3f left, right;
    std::tie(left, right) = tinynurbs::surfaceSplitV(srf, v);

    bool is_valid_l = tinynurbs::surfaceIsValid(left);
    bool is_valid_r = tinynurbs::surfaceIsValid(right);

    REQUIRE(is_valid_l == true);
    REQUIRE(is_valid_r == true);

    REQUIRE(left.degree_u == srf.degree_u);
    REQUIRE(right.degree_u == srf.degree_u);
    REQUIRE(left.degree_v == srf.degree_v);
    REQUIRE(right.degree_v == srf.degree_v);

    for (int i = 0; i < left.degree_v + 1; ++i) {
        int d = left.knots_v.size() - (left.degree_v + 1);
        REQUIRE(left.knots_v[d+i] == Approx(v));
    }

    for (int i = 0; i < right.degree_v + 1; ++i) {
        REQUIRE(right.knots_v[i] == Approx(v));
    }

    tinynurbs::surfaceSaveOBJ("left_rational_v.obj", left);
    tinynurbs::surfaceSaveOBJ("right_rational_v.obj", right);

    glm::vec3 pt1 = tinynurbs::surfacePoint(srf, left.knots_v[left.knots_v.size() - 1], 0.f);
    glm::vec3 pt2 = tinynurbs::surfacePoint(srf, right.knots_v[0], 0.f);
    REQUIRE(pt1.x == Approx(pt2.x));
    REQUIRE(pt1.y == Approx(pt2.y));
    REQUIRE(pt1.z == Approx(pt2.z));
}

TEST_CASE("surfaceSplitV (rational, interior split)", "[surface, rational, modify]")
{
    // The hemisphere is a Bezier in both directions (knots {0,0,0,0,1,1,1,1}),
    // but splitting at v = 0.5 (not a knot) forces full knot insertion of
    // multiplicity 3 before partitioning, exercising the non-degenerate path
    // of surfaceSplit along the v-direction (column partitioning).
    auto srf = getHemisphere();
    float v = 0.5f;
    tinynurbs::RationalSurface3f left, right;
    std::tie(left, right) = tinynurbs::surfaceSplitV(srf, v);

    REQUIRE(tinynurbs::surfaceIsValid(left) == true);
    REQUIRE(tinynurbs::surfaceIsValid(right) == true);
    REQUIRE(left.degree_u == srf.degree_u);
    REQUIRE(right.degree_u == srf.degree_u);
    REQUIRE(left.degree_v == srf.degree_v);
    REQUIRE(right.degree_v == srf.degree_v);

    // Knot-vector / control-point counts required by the validity relation.
    // Splitting at a non-knot parameter inserts it to full degree multiplicity
    // (r = degree), so the two pieces share one column: left + right = n + r + 1.
    REQUIRE(left.knots_v.size() == (left.degree_v + 1) + (size_t)left.control_points.cols());
    REQUIRE(right.knots_v.size() == (right.degree_v + 1) + (size_t)right.control_points.cols());
    REQUIRE((size_t)left.control_points.cols() + (size_t)right.control_points.cols() ==
            (size_t)srf.control_points.cols() + srf.degree_v + 1);
    REQUIRE((size_t)left.control_points.rows() == (size_t)srf.control_points.rows());
    REQUIRE((size_t)right.control_points.rows() == (size_t)srf.control_points.rows());

    // Left half ends at v with degree+1 copies, right half starts there
    for (int i = 0; i < left.degree_v + 1; ++i) {
        int d = left.knots_v.size() - (left.degree_v + 1);
        REQUIRE(left.knots_v[d + i] == Approx(v));
    }
    for (int i = 0; i < right.degree_v + 1; ++i) {
        REQUIRE(right.knots_v[i] == Approx(v));
    }

    // Geometric preservation along a u-isocurve for both halves
    const float u = 0.5f;
    const std::vector<float> sample_left = {0.f, 0.25f, v};
    for (float t : sample_left) {
        glm::vec3 p_orig = tinynurbs::surfacePoint(srf, u, t);
        glm::vec3 p_l = tinynurbs::surfacePoint(left, u, t);
        REQUIRE(p_l.x == Approx(p_orig.x));
        REQUIRE(p_l.y == Approx(p_orig.y));
        REQUIRE(p_l.z == Approx(p_orig.z));
    }
    const std::vector<float> sample_right = {v, 0.75f, 1.f};
    for (float t : sample_right) {
        glm::vec3 p_orig = tinynurbs::surfacePoint(srf, u, t);
        glm::vec3 p_r = tinynurbs::surfacePoint(right, u, t);
        REQUIRE(p_r.x == Approx(p_orig.x));
        REQUIRE(p_r.y == Approx(p_orig.y));
        REQUIRE(p_r.z == Approx(p_orig.z));
    }
}

TEST_CASE("surfaceReadOBJ and surfaceSaveOBJ (rational)", "[surface, obj]")
{
    auto srf = getHemisphere();
    
    tinynurbs::surfaceSaveOBJ("surface_rational.obj", srf);
    auto read_srf = tinynurbs::surfaceReadOBJ<float>("surface_rational.obj");
    
    REQUIRE(srf.degree_u == read_srf.degree_u);
    REQUIRE(srf.degree_v == read_srf.degree_v);
    REQUIRE(srf.knots_u.size() == read_srf.knots_u.size());
    for (int i = 0; i < srf.knots_u.size(); ++i) {
        REQUIRE(srf.knots_u[i] == Approx(read_srf.knots_u[i]));
    }
    REQUIRE(srf.knots_v.size() == read_srf.knots_v.size());
    for (int i = 0; i < srf.knots_v.size(); ++i) {
        REQUIRE(srf.knots_v[i] == Approx(read_srf.knots_v[i]));
    }
    REQUIRE(srf.control_points.rows() == read_srf.control_points.rows());
    REQUIRE(srf.control_points.cols() == read_srf.control_points.cols());
    for (int i = 0; i < srf.control_points.rows(); ++i) {
        for (int j = 0; j < srf.control_points.cols(); ++j) {
            REQUIRE(srf.control_points(i, j).x == Approx(read_srf.control_points(i, j).x));
            REQUIRE(srf.control_points(i, j).y == Approx(read_srf.control_points(i, j).y));
            REQUIRE(srf.control_points(i, j).z == Approx(read_srf.control_points(i, j).z));
        }
    }
    REQUIRE(srf.weights.rows() == read_srf.weights.rows());
    REQUIRE(srf.weights.cols() == read_srf.weights.cols());
    for (int i = 0; i < srf.weights.rows(); ++i) {
        for (int j = 0; j < srf.weights.cols(); ++j) {
            REQUIRE(srf.weights(i, j) == Approx(read_srf.weights(i, j)));
        }
    }
}