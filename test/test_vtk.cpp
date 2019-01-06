#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkQuad.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

#include <iostream>
#include <sstream>
using vtkPointsP                    = vtkSmartPointer<vtkPoints>;
using vtkUnstructuredGridP          = vtkSmartPointer<vtkUnstructuredGrid>;
using vtkQuadP                      = vtkSmartPointer<vtkQuad>;
using vtkXMLUnstructuredGridWriterP = vtkSmartPointer<vtkXMLUnstructuredGridWriter>;
using vtkIntArrayP                  = vtkSmartPointer<vtkIntArray>;
using Node    = std::vector<double>;
using Cell = std::vector<int>;                           

int main(int argc, char **argv)
{
    std::ostringstream fileName; 
    fileName << "vtk-test";
    auto writer = vtkXMLUnstructuredGridWriterP::New();
    fileName << "." << writer->GetDefaultFileExtension();
    writer->SetFileName((fileName.str()).c_str());
    auto dataSet = vtkUnstructuredGridP::New();
    auto pts = vtkPointsP::New();
    int num_pts = 9;
    pts->SetNumberOfPoints(num_pts);
    std::array<double, 27> points = {
        -1, -1, 0, 
         0, -1, 0,
         1, -1, 0,
        -1,  0, 0,
         0,  0, 0,
         1,  0, 0,
        -1,  1, 0,
         0,  1, 0,
         1,  1, 0
    };

    int id = 0;
    for(int i = 0; i < num_pts; i++)
    {
        pts->SetPoint(id, points[3*i + 0], points[3*i + 1], points[3*i + 2]);
        ++id;
    }

    std::vector<std::array<int, 4> >  cell = {
        {0, 1, 4, 3},
        {1, 2, 5, 4},
        {3, 4, 7, 6},
        {4, 5, 8, 7}
    };

    vtkQuadP quad = vtkQuadP::New();  // Assuming hex elements
    for(const auto & q : cell)
    {
        int num = 0;
        for(const auto & id : q)
        {
            (quad->GetPointIds())->SetId(num, id);
            ++num;
        }
        dataSet->InsertNextCell(quad->GetCellType(), quad->GetPointIds());
    } 

    std::vector<int> rank(4, 0);
    auto mpirank = vtkSmartPointer<vtkIntArray>::New();
    mpirank->SetName("mpirank");
    mpirank->SetArray(rank.data(), rank.size(), 1);
    dataSet->GetCellData()->AddArray(mpirank);

    dataSet->SetPoints(pts);
    // Remove unused memory
    dataSet->Squeeze();
    // Write the data
    writer->SetInputData(dataSet);
    writer->SetDataModeToAscii();
    writer->Write();
    return 0;
}

/*
void writeGrid(double time, const BoxArray& grid, std::ostringstream& fileName)
{
  // Create a writer
  auto writer = vtkXMLUnstructuredGridWriterP::New();
  // Append the default extension to the file name
  fileName << "." << writer->GetDefaultFileExtension();
  writer->SetFileName((fileName.str()).c_str());
  // Create a pointer to a VTK Unstructured Grid data set
  auto dataSet = vtkUnstructuredGridP::New();
  // Set up pointer to point data
  auto pts = vtkPointsP::New();
  // Count the total number of points to be saved
  int num_pts = grid.getNumberOfPoints(); // Implementation is user-dependent
  pts->SetNumberOfPoints(num_pts);
  // Add the time
  addTimeToVTKDataSet(time, dataSet);
  // Get the nodes and elements that are used to describe the grid
  std::vector<Node>    nodes = grid.getNodes();       // Implementation is user-dependent
  std::vector<Element> elements = grid.getElements(); // Implementation is user-dependent
  // Add the processor boundaries to the unstructured grid cell data
  addElementsToVTKDataSet(nodes, elements, pts, dataSet);
  // Set the points
  dataSet->SetPoints(pts);
  // Remove unused memory
  dataSet->Squeeze();
  // Write the data
  writer->SetInput(dataSet);
  writer->SetDataModeToAscii();
  writer->Write();
}

OutputVTK::addTimeToVTKDataSet(double time,
                               vtkUnstructuredGridP& dataSet)
{
  auto array = vtkDoubleArrayP::New();
  array->SetName("TIME");
  array->SetNumberOfTuples(1);
  array->SetTuple1(0, time);
  dataSet->GetFieldData()->AddArray(array);
}


void
OutputVTK::addElementsToVTKDataSet(const std::vector<Node>& nodes,
                                   const std::vector<Element>& elements,
                                   vtkPointsP& pts,
                                   vtkUnstructuredGridP& dataSet)
{
  // Set the coordinates of the nodes
  int id = 0;
  for (const auto& node : nodes) {
    pts->SetPoint(id, node[0], node[1], node[2]);
    ++id;
  }
  // Set the element connectivities
  auto hex = vtkHexahedronP::New();  // Assuming hex elements
  for (const auto& element : elements) {
    // Get node ids and assign them to a hex element (ids start from 1)
    int nodeNum = 0;
    for (const auto& id : element) {
      hex->GetPointIds()->SetId(nodeNum, id - 1);
      ++nodeNum;
    }
    dataSet->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
  }
}*/
