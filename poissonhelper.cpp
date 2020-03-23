#include "poissonhelper.h"
#include <iostream>
#include <fstream>
#include <string>

#include <libigl/include/igl/copyleft/marching_cubes.h>

PoissonHelper::PoissonHelper()
{
}

int PoissonHelper::cubeToRowIndex(int x, int y, int z, int gridWidth, int gridHeight){
    return x + gridWidth * y + gridWidth * gridHeight * z;
}

// Approximation to the inner product of <F, laplacian(F')>
float PoissonHelper::integral_f_dd_fPrime(float normalizedRadius)
{
    if(normalizedRadius >= -1 && normalizedRadius < 0){
        return 0.0f;
    }
    if(normalizedRadius >= 0.0f && normalizedRadius < 1.0f){
        return -4.0f / 3.0f;
    }
    if(normalizedRadius >= 1.0f && normalizedRadius < 2){
        return 2.0f / 3.0f;
    }
    return 0.0f;
}

// Approximation to the inner product of <F, gradient(F')>
float PoissonHelper::integral_f_d_fPrime(float normalizedRadius)
{
    if(normalizedRadius >= -1 && normalizedRadius < 0){
        return 0.0f;
    }
    if(normalizedRadius >= 0.0f && normalizedRadius < 1.0f){
        return 1.0f;
    }
    if(normalizedRadius >= 1.0f && normalizedRadius < 2){
        return -1.0f;
    }
    return 0.0f;
}

// Approximation to the inner product of <F, F'>
float PoissonHelper::integral_f_fPrime(float normalizedRadius)
{
    if(normalizedRadius >= -1 && normalizedRadius < 0){
        return 0.0;
    }
    if(normalizedRadius >= 0.0f && normalizedRadius < 1.0f){
        return 0.45f;
    }
    if(normalizedRadius >= 1.0f && normalizedRadius < 2){
        return 0.775f;
    }
    return 0.0f;
}

// The basis method we'll be using. Feel free to experiment with other
// bases.
float PoissonHelper::basis(float normalizedRadius)
{
    if(normalizedRadius >= -1 && normalizedRadius < 0){
        return 0.5f + normalizedRadius + 0.5f * normalizedRadius * normalizedRadius;
    }
    if(normalizedRadius >= 0.0f && normalizedRadius < 1.0f){
        return 0.5f + normalizedRadius - normalizedRadius * normalizedRadius;
    }
    if(normalizedRadius >= 1.0f && normalizedRadius < 2){
        return 2.0f - 2.0f * normalizedRadius + 0.5f * normalizedRadius * normalizedRadius;
    }
    return 0.0f;
}

// Fill look-up tables with inner products.
void PoissonHelper::fillLookupTables(int gridWidth, int gridHeight, int gridDepth)
{
    for(int i = 0; i < gridWidth; i++){
        for (int j = 0; j < gridHeight; j++){
            for (int k = 0; k < gridDepth; k ++ ){
                int index = cubeToRowIndex(i, j, k, gridWidth, gridHeight);

                for(int m = -1; m < 2; m++){
                    for(int n = -1; n < 2; n++){
                        for(int z = -1;z < 2; z++){
                            int neighbor = cubeToRowIndex(i + m, j + n, z + k, gridWidth, gridHeight);
                            if(neighbor < 0 || i + m >= gridWidth || j + n >= gridHeight ||z + k >= gridDepth){
                                 continue;
                            }
                            f_dot_fprime_x[index][neighbor] += integral_f_fPrime(m);
                            f_dot_dd_fprime_x[index][neighbor] += integral_f_dd_fPrime(m);
                            f_dot_d_fprime_x[index][neighbor] += integral_f_d_fPrime( m);

                            f_dot_fprime_y[index][neighbor] += integral_f_fPrime(n);
                            f_dot_dd_fprime_y[index][neighbor] += integral_f_dd_fPrime(n);
                            f_dot_d_fprime_y[index][neighbor] += integral_f_d_fPrime( n);

                            f_dot_fprime_z[index][neighbor] += integral_f_fPrime(z);
                            f_dot_dd_fprime_z[index][neighbor] += integral_f_dd_fPrime(z);
                            f_dot_d_fprime_z[index][neighbor] += integral_f_d_fPrime( z);
                        }
                    }
                }
            }
        }
    }
}

// Use libigl to mesh the implicit function, given an isovalue
void PoissonHelper::marchingCubes(Eigen::VectorXf surface, float isoValue, int voxelGridWidth, int voxelGridHeight,
                                  int voxelGridDepth, float cellWidth, Vector3f corner, Eigen::MatrixXf &V, Eigen::MatrixXi &F)
{

    Eigen::MatrixXf points(voxelGridDepth * voxelGridHeight * voxelGridWidth, 3);

    for(int i = 0; i < voxelGridWidth; i++){
        for(int j = 0; j < voxelGridHeight; j++){
            for(int k = 0; k < voxelGridDepth; k++){
                points(cubeToRowIndex(i, j, k, voxelGridWidth, voxelGridHeight), 0) = corner[0] + cellWidth * i;
                points(cubeToRowIndex(i, j, k, voxelGridWidth, voxelGridHeight), 1) = corner[1] + cellWidth * j;
                points(cubeToRowIndex(i, j, k, voxelGridWidth, voxelGridHeight), 2) = corner[2] + cellWidth * k;
            }
        }
    }

    igl::copyleft::marching_cubes(surface, points, voxelGridWidth, voxelGridHeight, voxelGridDepth, isoValue, V, F);
}

// load in a point cloud
void PoissonHelper::loadPointsFromFile(const std::string &filePath, std::vector<Vector3f> &positions, std::vector<Vector3f> &normals)
{
    std::ifstream pointFile;
    pointFile.open(filePath.c_str());

    if(pointFile.is_open()){
        string data;
        bool past_header = false;

        while (getline(pointFile, data)) {
            if(!past_header && strcmp(data.c_str(), "end_header")){
                past_header = true;
            }
            if(!past_header){
                continue;
            }

            std::vector<string> lineData;

            char *dup = strdup(data.c_str());
            char * split = strtok(dup, " ");
            while(split != NULL)
            {
                lineData.push_back(split);
                split = strtok(NULL, " ");
            }
            if(lineData.size() < 6){
                continue;
            }
            positions.push_back( Vector3f(stof(lineData[0]), stof(lineData[1]), stof(lineData[2])));
            normals.push_back( Vector3f(stof(lineData[3]), stof(lineData[4]), stof(lineData[5])));
        }
        std::cout << "Loaded data" << std::endl;
    } else {
        std::cerr << "Failed to load/parse input file" << std::endl;
        return;
    }
}

// You get this for free! Given a point cloud, determine the dimensions
// of the uniform grid, the size of each individual cell, and the
// lower/left/back corner of the entire uniform grid.
void PoissonHelper::getGridDimensions(const std::vector<Vector3f> &positions, int &voxelGridWidth, int &voxelGridHeight,
                                      int &voxelGridDepth, float &cellWidth, Vector3f &corner)
{
    float maxX = -std::numeric_limits<float>::infinity();
    float maxY = -std::numeric_limits<float>::infinity();
    float maxZ = -std::numeric_limits<float>::infinity();

    float minX = std::numeric_limits<float>::infinity();
    float minY = std::numeric_limits<float>::infinity();
    float minZ = std::numeric_limits<float>::infinity();

    for(int i = 0; i < positions.size(); i++){
        if(positions[i][0] > maxX){
            maxX = positions[i][0];
        }

        if(positions[i][1] > maxY){
            maxY = positions[i][1];
        }

        if(positions[i][2] > maxZ){
            maxZ = positions[i][2];
        }

        if(positions[i][0] < minX){
            minX = positions[i][0];
        }

        if(positions[i][1] < minY){
            minY = positions[i][1];
        }

        if(positions[i][2] < minZ){
            minZ = positions[i][2];
        }
    }
    std::cout << "Bounding Box: " << minX << " " << maxX << " " << minY << " " << maxY << " " << minZ << " " << maxZ <<" " << positions.size() << std::endl;
    // number of extra grid cells on each side of the voxel cube
    float padding = 8.0;
    // 20 is an arbitrary number that can be thought of as the granularity of the
    // uniform grid. As this number increases, the grid will become denser, and the
    // program will take longer to run. Unless, of course, you'd like to try an octree...
    cellWidth = (maxX - minX) / 20.0f;

    voxelGridWidth = ceil( ((maxX - minX) + 2 * padding * cellWidth) / cellWidth);
    voxelGridHeight = ceil( ((maxY - minY) + 2 * padding * cellWidth) / cellWidth);
    voxelGridDepth = ceil( ((maxZ - minZ) + 2 * padding * cellWidth) / cellWidth) ;

    corner = Vector3f(minX - padding * cellWidth, minY - padding * cellWidth, minZ - padding * cellWidth);
}

// You get this for free! Given vertices and faces, write the mesh out to file.
void PoissonHelper::saveAsMesh(const std::string &filePath,const MatrixXf &vertices,
                      const MatrixXi &faces)
{
    std::ofstream outfile;
    outfile.open(filePath);

    // Write vertices
    for (size_t i = 0; i < vertices.rows(); i++)
    {
        outfile << "v " << vertices(i,0) << " " << vertices(i,1) << " " << vertices(i,2) << std::endl;
    }

    // Write faces
    for (size_t i = 0; i < faces.rows(); i++)
    {
        outfile << "f " << (faces(i,0) + 1) << " " << (faces(i, 1) + 1) << " " << (faces(i,2) + 1) << std::endl;
    }

    outfile.close();
}
