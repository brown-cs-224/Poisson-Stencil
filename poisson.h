#ifndef POISSON_H
#define POISSON_H

#include <vector>

#include <Eigen/StdVector>
#include <Eigen/Sparse>
#include "poissonhelper.h"
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix2f)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3f)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3i)

using namespace Eigen;
class Poisson
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Poisson();
    void reconstruct(const std::string infile, const std::string outfile);

    void fillDivergenceVector(std::vector<Vector3f> positions, std::vector<Vector3f> normals,
                              int gridWidth, int gridHeight, int gridDepth, float cellWidth,
                              Vector3f corner, std::map<int, Vector3f> vectorField, Eigen::VectorXf &b);

    void pointCloudToVectorField(std::vector<Vector3f> positions, std::vector<Vector3f> normals,
                                          int gridWidth, int gridHeight, int gridDepth, float cellWidth,
                                          Vector3f corner, std::map<int, Vector3f> &vectorField);

    void fillLaplacianMatrix(int gridWidth, int gridHeight, int gridDepth, Eigen::SparseMatrix<float> &L);

    void meshIsosurface(Eigen::VectorXf surface, std::vector<Vector3f> positions,
                        int gridWidth, int gridHeight, int gridDepth, float cellWidth,
                        Vector3f corner, std::string outfile);

private:
    PoissonHelper poissonHelper;
};

#endif // POISSON_H
