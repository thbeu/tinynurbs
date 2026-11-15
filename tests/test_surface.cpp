#include <tinynurbs/tinynurbs.h>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <cmath>
#include "catch.hpp"

using namespace std;

tinynurbs::Surface3f getBilinearPatch() {
    tinynurbs::Surface3f srf;
    srf.degree_u = 1;
    srf.degree_v = 1;
    srf.knots_u = {0, 0, 1, 1};
    srf.knots_v = {0, 0, 1, 1};
    // 2x2 grid (tinynurbs::array2) of control points
    srf.control_points = {2, 2, 
                          {glm::vec3(-1, 0, 1), glm::vec3(-1, 0, -1),
                           glm::vec3(1, 0, 1), glm::vec3(1, 0, -1)
                          }
                         };
    return srf;
}

TEST_CASE("surfacePoint (non-rational)", "[surface, non-rational, evaluate]")
{
    auto srf = getBilinearPatch();
    glm::vec3 pt1 = tinynurbs::surfacePoint(srf, 0.f, 0.f);
    glm::vec3 pt2 = tinynurbs::surfacePoint(srf, 0.5f, 0.5f);
    REQUIRE(pt1.x == Approx(-1));
    REQUIRE(pt1.y == Approx(0));
    REQUIRE(pt1.z == Approx(1));
    REQUIRE(pt2.x == Approx(0));
    REQUIRE(pt2.y == Approx(0));
    REQUIRE(pt2.z == Approx(0));
}

TEST_CASE("surfaceTangent (non-rational)", "[surface, non-rational, evaluate]")
{
    auto srf = getBilinearPatch();
    glm::vec3 tgt_u, tgt_v;
    std::tie(tgt_u, tgt_v) = tinynurbs::surfaceTangent(srf, 0.5f, 0.25f);
    REQUIRE(glm::length(tgt_u) == Approx(1));
    REQUIRE(tgt_u.x == Approx(1));
    REQUIRE(tgt_u.y == Approx(0));
    REQUIRE(tgt_u.z == Approx(0));
    REQUIRE(glm::length(tgt_v) == Approx(1));
    REQUIRE(tgt_v.x == Approx(0));
    REQUIRE(tgt_v.y == Approx(0));
    REQUIRE(tgt_v.z == Approx(-1));
}


TEST_CASE("surfaceNormal (non-rational)", "[surface, non-rational, evaluate]")
{
    auto srf = getBilinearPatch();
    glm::vec3 n = tinynurbs::surfaceNormal(srf, 0.5f, 0.5f);
    REQUIRE(glm::length(n) == Approx(1));
    REQUIRE(n.x == Approx(0));
    REQUIRE(n.y == Approx(-1));
    REQUIRE(n.z == Approx(0));

    n = tinynurbs::surfaceNormal(srf, 0.25f, 0.75f);
    REQUIRE(glm::length(n) == Approx(1));
    REQUIRE(n.x == Approx(0));
    REQUIRE(n.y == Approx(-1));
    REQUIRE(n.z == Approx(0));
}

TEST_CASE("surfaceIsValid (non-rational)", "[surface, check]")
{
    auto srf = getBilinearPatch();
    bool is_valid = tinynurbs::surfaceIsValid(srf);
    REQUIRE(is_valid == true);

    srf = getBilinearPatch();
    srf.degree_u = 5;
    is_valid = tinynurbs::surfaceIsValid(srf);
    REQUIRE(is_valid == false);

    srf = getBilinearPatch();
    srf.degree_v = 4;
    is_valid = tinynurbs::surfaceIsValid(srf);
    REQUIRE(is_valid == false);
}

TEST_CASE("surfaceInsertKnotU (non-rational)", "[surface, non-rational, modify]")
{
    auto srf = getBilinearPatch();
    unsigned int repeat = 1;

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

TEST_CASE("surfaceInsertKnotV (non-rational)", "[surface, non-rational, modify]")
{
    auto srf = getBilinearPatch();
    unsigned int repeat = 1;

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

TEST_CASE("surfaceSplitU (non-rational)", "[surface, non-rational, modify]")
{
    auto srf = getBilinearPatch();
    float u = 0.25f;
    tinynurbs::Surface3f left, right;
    std::tie(left, right) = tinynurbs::surfaceSplitU(srf, u);

    bool is_valid_l = tinynurbs::surfaceIsValid(left);
    bool is_valid_r = tinynurbs::surfaceIsValid(right);

    REQUIRE(is_valid_l == true);
    REQUIRE(is_valid_r == true);

    REQUIRE(left.degree_u == srf.degree_u);
    REQUIRE(right.degree_u == srf.degree_u);
    REQUIRE(left.degree_v == srf.degree_v);
    REQUIRE(right.degree_v == srf.degree_v);

    for (unsigned int i = 0; i < left.degree_u + 1; ++i) {
        size_t d = left.knots_u.size() - (left.degree_u + 1);
        REQUIRE(left.knots_u[d+i] == Approx(u));
    }

    for (unsigned int i = 0; i < right.degree_u + 1; ++i) {
        REQUIRE(right.knots_u[i] == Approx(u));
    }

    tinynurbs::surfaceSaveOBJ("left_nonrational_u.obj", left);
    tinynurbs::surfaceSaveOBJ("right_nonrational_u.obj", right);

    glm::vec3 pt1 = tinynurbs::surfacePoint(srf, left.knots_u[left.knots_u.size() - 1], 0.f);
    glm::vec3 pt2 = tinynurbs::surfacePoint(srf, right.knots_u[0], 0.f);
    REQUIRE(pt1.x == Approx(pt2.x));
    REQUIRE(pt1.y == Approx(pt2.y));
    REQUIRE(pt1.z == Approx(pt2.z));
}

TEST_CASE("surfaceSplitV (non-rational)", "[surface, non-rational, modify]")
{
    auto srf = getBilinearPatch();
    float v = 0.25f;
    tinynurbs::Surface3f left, right;
    std::tie(left, right) = tinynurbs::surfaceSplitV(srf, v);

    bool is_valid_l = tinynurbs::surfaceIsValid(left);
    bool is_valid_r = tinynurbs::surfaceIsValid(right);

    REQUIRE(is_valid_l == true);
    REQUIRE(is_valid_r == true);

    REQUIRE(left.degree_u == srf.degree_u);
    REQUIRE(right.degree_u == srf.degree_u);
    REQUIRE(left.degree_v == srf.degree_v);
    REQUIRE(right.degree_v == srf.degree_v);

    for (unsigned int i = 0; i < left.degree_v + 1; ++i) {
        size_t d = left.knots_v.size() - (left.degree_v + 1);
        REQUIRE(left.knots_v[d+i] == Approx(v));
    }

    for (unsigned int i = 0; i < right.degree_v + 1; ++i) {
        REQUIRE(right.knots_v[i] == Approx(v));
    }

    tinynurbs::surfaceSaveOBJ("left_nonrational_v.obj", left);
    tinynurbs::surfaceSaveOBJ("right_nonrational_v.obj", right);

    glm::vec3 pt1 = tinynurbs::surfacePoint(srf, left.knots_v[left.knots_v.size() - 1], 0.f);
    glm::vec3 pt2 = tinynurbs::surfacePoint(srf, right.knots_v[0], 0.f);
    REQUIRE(pt1.x == Approx(pt2.x));
    REQUIRE(pt1.y == Approx(pt2.y));
    REQUIRE(pt1.z == Approx(pt2.z));
}

TEST_CASE("surfaceSplitU (non-rational, non-bezier)", "[surface, non-rational, modify]")
{
    // Surface with a cubic (degree 3) u-direction having one interior knot
    // and a linear v-direction. 5x2 grid of control points.
    tinynurbs::Surface3f srf;
    srf.degree_u = 3;
    srf.degree_v = 1;
    srf.knots_u = {0, 0, 0, 0, 0.5f, 1, 1, 1, 1};
    srf.knots_v = {0, 0, 1, 1};
    srf.control_points = {5, 2,
                          {glm::vec3(0, 0, 0),  glm::vec3(0, 0, -1),
                           glm::vec3(1, 1, 0),  glm::vec3(1, 1, -1),
                           glm::vec3(2, 0, 0),  glm::vec3(2, 0, -1),
                           glm::vec3(3, 0, 0),  glm::vec3(3, 0, -1),
                           glm::vec3(4, 0, 0),  glm::vec3(4, 0, -1)
                          }
                         };

    // Split at u = 0.25, which is NOT an existing knot (forces full knot
    // insertion before partitioning) in the cubic u-direction.
    float u = 0.25f;
    tinynurbs::Surface3f left, right;
    std::tie(left, right) = tinynurbs::surfaceSplitU(srf, u);

    REQUIRE(tinynurbs::surfaceIsValid(left) == true);
    REQUIRE(tinynurbs::surfaceIsValid(right) == true);
    REQUIRE(left.degree_u == srf.degree_u);
    REQUIRE(right.degree_u == srf.degree_u);
    REQUIRE(left.degree_v == srf.degree_v);
    REQUIRE(right.degree_v == srf.degree_v);

    // Knot-vector / control-point counts required by the validity relation.
    // Splitting at a non-knot parameter inserts it to full degree multiplicity
    // (r = degree), so the two pieces share one row: left + right = n + r + 1.
    REQUIRE(left.knots_u.size() == (left.degree_u + 1) + (size_t)left.control_points.rows());
    REQUIRE(right.knots_u.size() == (right.degree_u + 1) + (size_t)right.control_points.rows());
    REQUIRE((size_t)left.control_points.rows() + (size_t)right.control_points.rows() ==
            (size_t)srf.control_points.rows() + srf.degree_u + 1);
    REQUIRE((size_t)left.control_points.cols() == (size_t)srf.control_points.cols());
    REQUIRE((size_t)right.control_points.cols() == (size_t)srf.control_points.cols());

    // Left half ends at u with degree+1 copies, right half starts there
    for (unsigned int i = 0; i < left.degree_u + 1; ++i) {
        size_t d = left.knots_u.size() - (left.degree_u + 1);
        REQUIRE(left.knots_u[d + i] == Approx(u));
    }
    for (unsigned int i = 0; i < right.degree_u + 1; ++i) {
        REQUIRE(right.knots_u[i] == Approx(u));
    }

    // Geometric preservation along a v-isocurve for both halves
    const float v = 0.5f;
    const std::vector<float> sample_left = {0.f, 0.1f, u};
    for (float t : sample_left) {
        glm::vec3 p_orig = tinynurbs::surfacePoint(srf, t, v);
        glm::vec3 p_l = tinynurbs::surfacePoint(left, t, v);
        REQUIRE(p_l.x == Approx(p_orig.x));
        REQUIRE(p_l.y == Approx(p_orig.y));
        REQUIRE(p_l.z == Approx(p_orig.z));
    }
    const std::vector<float> sample_right = {u, 0.75f, 1.f};
    for (float t : sample_right) {
        glm::vec3 p_orig = tinynurbs::surfacePoint(srf, t, v);
        glm::vec3 p_r = tinynurbs::surfacePoint(right, t, v);
        REQUIRE(p_r.x == Approx(p_orig.x));
        REQUIRE(p_r.y == Approx(p_orig.y));
        REQUIRE(p_r.z == Approx(p_orig.z));
    }
}

TEST_CASE("surfaceReadOBJ and surfaceSaveOBJ (non-rational)", "[surface, non-rational, obj]")
{
    auto srf = getBilinearPatch();
    
    tinynurbs::surfaceSaveOBJ("surface_nonrational.obj", srf);
    auto read_srf = tinynurbs::surfaceReadOBJ<float>("surface_nonrational.obj");
    
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
}