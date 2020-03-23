#ifndef POISSONHELPER_H
#define POISSONHELPER_H

#include <vector>

#include <Eigen/StdVector>
#include <Eigen/Sparse>
using namespace Eigen;
using namespace std;

class PoissonHelper
{
public:
    PoissonHelper();
    int cubeToRowIndex(int x, int y, int z, int gridWidth, int gridHeight);

    float basis(float normalizedRadius);
    void fillLookupTables(int gridWidth, int gridHeight, int gridDepth);
    void loadPointsFromFile(const std::string &filePath, std::vector<Eigen::Vector3f> &positions, std::vector<Eigen::Vector3f> &normals);
    void getGridDimensions(const std::vector<Vector3f> &positions, int &voxelGridWidth, int &voxelGridHeight,
                           int &voxelGridDepth, float &cellWidth, Vector3f &corner);

    void marchingCubes(Eigen::VectorXf surface, float isoValue, int voxelGridWidth, int voxelGridHeight,
                       int voxelGridDepth, float cellWidth, Vector3f corner, Eigen::MatrixXf &V, Eigen::MatrixXi &F);

    void saveAsMesh(const std::string &filePath,const MatrixXf &vertices,
                          const MatrixXi &faces);

    float integral_f_dd_fPrime(float normalizedRadius);
    float integral_f_d_fPrime(float normalizedRadius);
    float integral_f_fPrime(float normalizedRadius);

    std::map<int, std::map<int, float>> f_dot_fprime_x;
    std::map<int, std::map<int, float>> f_dot_d_fprime_x;
    std::map<int, std::map<int, float>> f_dot_dd_fprime_x;

    std::map<int, std::map<int, float>> f_dot_fprime_y;
    std::map<int, std::map<int, float>> f_dot_d_fprime_y;
    std::map<int, std::map<int, float>> f_dot_dd_fprime_y;

    std::map<int, std::map<int, float>> f_dot_fprime_z;
    std::map<int, std::map<int, float>> f_dot_d_fprime_z;
    std::map<int, std::map<int, float>> f_dot_dd_fprime_z;
};

#endif // POISSONHELPER_H
