#include "poisson.h"

#include <iostream>
#include <fstream>

#include <QFileInfo>
#include <QString>

#include <libigl/include/igl/copyleft/marching_cubes.h>

#define TINYOBJLOADER_IMPLEMENTATION
#include "util/tiny_obj_loader.h"
#include <Eigen/IterativeLinearSolvers>
using namespace Eigen;
using namespace std;

Poisson::Poisson()
{
}

void Poisson::reconstruct(const std::string infile, const std::string outfile)
{
    std::vector<Vector3f> positions;
    std::vector<Vector3f> normals;
    // Load in point cloud to fill an ordered list of positions and normals
    poissonHelper.loadPointsFromFile(infile, positions, normals);

    float cellWidth = 0.0f;
    int gridWidth = 0;
    int gridHeight = 0;
    int gridDepth = 0;
    Vector3f corner = Vector3f(0.0,0.0,0.0);
    // Populate variables for the dimensions of the grid, width of each cell, and lower-left-back corner of the grid.
    poissonHelper.getGridDimensions(positions, gridWidth, gridHeight, gridDepth, cellWidth, corner);

    // Fill inner-product look up tables.
    poissonHelper.fillLookupTables(gridWidth, gridHeight, gridDepth);

    // Use the helper methods below to fill these variables.
    Eigen::SparseMatrix<float> L;
    Eigen::VectorXf b;
    std::map<int, Vector3f> vectorField;
    fillLaplacianMatrix(gridWidth, gridHeight, gridDepth, L);
    pointCloudToVectorField(positions, normals, gridWidth, gridHeight, gridDepth,cellWidth, corner, vectorField);
    fillDivergenceVector(positions, normals, gridWidth, gridHeight, gridDepth, cellWidth, corner, vectorField, b);

    // Use Eigen's sparse library to solve Lx = b!
    Eigen::VectorXf x;

    // Calculate an isovalue and mesh the surface.
    meshIsosurface(x, positions, gridWidth, gridHeight, gridDepth, cellWidth, corner, outfile);

}

/**
 * @brief Fill the Eigen Sparse Matrix L, the Laplacian.
 * @param gridWidth, gridHeight, gridDepth: dimensions of the uniform grid.
 * @param L: fill this in!
 */
void Poisson::fillLaplacianMatrix(int gridWidth, int gridHeight, int gridDepth, Eigen::SparseMatrix<float> &L)
{

    // If the uniform grid is gridWidth nodes wide, gridHeight nodes tall, and gridDepth
    // nodes deep, the total number of nodes is...?
    int numberOfNodes = 0; //replace me!
    L.resize(numberOfNodes, numberOfNodes);
    L.setZero();

    // Fill the matrix with the Laplacian of the vector field. You can retrieve the
    // necessary inner products for each dimension from the look up tables:
    // x dimension: poissonHelper.f_dot_fprime_x, poissonHelper.f_dot_dd_fprime_x
    // y dimension: poissonHelper.f_dot_fprime_y, poissonHelper.f_dot_dd_fprime_y
    // z dimension: poissonHelper.f_dot_fprime_z, poissonHelper.f_dot_dd_fprime_z

    // Example:
    // poissonHelper.f_dot_fprime_x[1D index of node1][1D index of node2] (inner product of function 1 with function 2)
    // poissonHelper.f_dot_dd_fprime_x[1D index of node1][1D index of node2] (inner product of function 1 with 2nd derivative of function 2)

    // Check the handout for how to combine these inner products to find the value at each entry of L! Also, remember
    // that most of these values will be 0, unless node1 and node2 are within 2 neighbors of each other...
}

/**
 * @brief Given an oriented point cloud, convert the point cloud into a vector field discretized
 * over a uniform grid.
 * @param positions: a list of point positions
 * @param normals: a list of normals for every point
 * @param gridWidth, gridHeight, gridDepth: dimensions of the uniform grid
 * @param cellWidth: size of each individual cell in the grid
 * @param corner: the lower/left/back corner of the entire uniform grid. This will come in
 * handy for converting positions to indices in the uniform grid.
 * @param vectorField: fill this in!
 */
void Poisson::pointCloudToVectorField(std::vector<Vector3f> positions, std::vector<Vector3f> normals,
                                      int gridWidth, int gridHeight, int gridDepth, float cellWidth,
                                      Vector3f corner, std::map<int, Vector3f> &vectorField)
{

    // If you're on step 1: Ignore this method.

    // If you're on step 2: fill out the vectorField map in the following format:
    // vectorField[1D index of node1] = weighted sum of normals within node1's neighborhood

    // Iterate through all the points, and distribute every point's normal to its 8 surrounding
    // cells. The weights are determined by the neighbor's basis functions, given a point's distance
    // from the neighbor's cell center.
    // You can get the weight with:
    // weight = poissonHelper.basis(normalized_x_distance_from_node) * poissonHelper.basis(normalized_y_distance_from_node) * poissonHelper.basis(normalized_z_distance_from_node)


    for(int i = 0; i < positions.size(); i++){
        // Here's a start: calculate the grid indices of each point from the position.
        int xIndex = floor( (positions[i][0] - corner[0]) / cellWidth);
        int yIndex = floor( (positions[i][1] - corner[1]) / cellWidth);
        int zIndex = floor( (positions[i][2] - corner[2]) / cellWidth);
        Vector3f grid_position = positions[i] - corner;
    }
}

/**
 * @brief Fill out the divergence vector (Step 2 only).
 * @param positions: a list of point positions
 * @param normals: a list of normals for every point
 * @param gridWidth, gridHeight, gridDepth: dimension of the uniform grid
 * @param cellWidth: size of each individual cell in the grid
 * @param corner: the lower/left/back corner of the entire uniform grid. This will come in
 * handy for converting positions to indices in the uniform grid.
 * @param vectorField: result of pointCloudToVectorField()
 * @param b: fill this out!
 */
void Poisson::fillDivergenceVector(std::vector<Vector3f> positions, std::vector<Vector3f> normals,
                                   int gridWidth, int gridHeight, int gridDepth, float cellWidth,
                                   Vector3f corner, std::map<int, Vector3f> vectorField, Eigen::VectorXf &b)
{
    int numNodes = gridWidth * gridHeight * gridDepth;
    b.resize(numNodes);
    b.setZero();

    // if you're on step 1, IGNORE THIS METHOD. Otherwise, you'll need to replace this for loop.
    for(int i = 0; i < positions.size(); i++){
        // Shift the positions into a positive coordinate frame, then convert to indices on the
        // uniform grid.
        int index_x = floor( (positions[i][0] - corner[0]) / cellWidth);
        int index_y = floor( (positions[i][1] - corner[1]) / cellWidth);
        int index_z = floor( (positions[i][2] - corner[2]) / cellWidth);
        // find the associated flattened index of the cell
        int index = poissonHelper.cubeToRowIndex(index_x, index_y, index_z, gridWidth, gridHeight);
        // set a constant divergence
        b[index] = 1.0f;
    }


    // Step 2:
    // Remove the above for loop.
    // Fill the vector b with the Divergence of the vector field. You can retrieve the
    // necessary inner products for each dimension from the look up tables:
    // x dimension: poissonHelper.f_dot_fprime_x, poissonHelper.f_dot_d_fprime_x
    // y dimension: poissonHelper.f_dot_fprime_y, poissonHelper.f_dot_d_fprime_y
    // z dimension: poissonHelper.f_dot_fprime_z, poissonHelper.f_dot_d_fprime_z

    // Example:
    // poissonHelper.f_dot_fprime_x[1D index of node1][1D index of node2] (inner product of function 1 with function 2)
    // poissonHelper.f_dot_d_fprime_x[1D index of node1][1D index of node2] (inner product of function 1 with 1st derivative of function 2)

    // Check the handout for how to combine this information to find the gradient of the basis functions!
    // Remember that you will need to take the dot product of the gradient with the vector field.

}

/**
 * @brief Calculate the isovalue and mesh the smoothed implicit function
 * @param weights: the x vector that was the result of solving Lx = b
 * @param positions: a list of point positions
 * @param normals: a list of normals for every point
 * @param gridWidth, gridHeight, gridDepth: dimension of the uniform grid
 * @param cellWidth: size of each individual cell in the grid
 * @param corner: the lower/left/back corner of the entire uniform grid. This will come in
 * handy for converting positions to indices in the uniform grid.
 * @param outfile: output .obj file
 */
void Poisson::meshIsosurface(Eigen::VectorXf weights, std::vector<Vector3f> positions,
                             int gridWidth, int gridHeight, int gridDepth, float cellWidth,
                             Vector3f corner, std::string outfile)
{
    // Calculate an appropriate isovalue!
    float isoValue = 0.0f;

    for(int i = 0; i < positions.size(); i++){
        // Here's a start: calculate the grid indices of each point from the position.
        int xIndex = floor( (positions[i][0] - corner[0]) / cellWidth);
        int yIndex = floor( (positions[i][1] - corner[1]) / cellWidth);
        int zIndex = floor( (positions[i][2] - corner[2]) / cellWidth);
        Vector3f grid_position = positions[i] - corner;
    }


    std::cout << "Using isovalue of " << isoValue << std::endl;

    // The stencil will take it from here. Given a good isovalue and the weights vector, we
    // execute libigl's marching cubes method and write out the vertices and faces.
    if(weights.size() > 0){
        Eigen::MatrixXf vertices;
        Eigen::MatrixXi faces;
        poissonHelper.marchingCubes(weights, isoValue, gridWidth, gridHeight, gridDepth, cellWidth, corner, vertices, faces);
        poissonHelper.saveAsMesh(outfile, vertices, faces);
    } else {
        std::cout << "No surface provided to marching cubes!" << std::endl;
    }
}
