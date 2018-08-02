#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkVersion.h>
#include "collid.h"

void vtkplotVectorSurface(std::vector<CD_HSE*>& hse_list, const char* fname)
{
  vtkSmartPointer<vtkPoints> pts =
    vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkCellArray> cells =
    vtkSmartPointer<vtkCellArray>::New();

  vtkSmartPointer<vtkUnsignedCharArray> colors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("CollsnRegion");

  vtkSmartPointer<vtkFloatArray> impulses =
    vtkSmartPointer<vtkFloatArray>::New();
  impulses->SetNumberOfComponents(3);
  impulses->SetName("CollsnImpulse");

  unsigned char red[3] = {255, 0, 0};
  unsigned char green[3] = {0, 255, 0};

  int id = 0;
  unsortHseList(hse_list);
  for (unsigned i = 0; i < hse_list.size(); ++i)
  for (int j = 0; j < hse_list[i]->num_pts(); ++j){
	POINT* p = hse_list[i]->Point_of_hse(j);
	if (sorted(p)) continue;
	p->indx = id++;
	sorted(p) = YES;
  }
  for (unsigned i = 0; i < hse_list.size(); ++i)
  {
         vtkSmartPointer<vtkPolygon>cell =
         vtkSmartPointer<vtkPolygon>::New();
	 cell->GetPointIds()->SetNumberOfIds(hse_list[i]->num_pts());
	 
	 bool has_collsn = false;
     	 for (int j = 0; j < hse_list[i]->num_pts(); ++j){
	     //insert geometry
	     POINT* p = hse_list[i]->Point_of_hse(j);
	     id = p->indx;
             pts->InsertPoint(id,Coords(p));
	     cell->GetPointIds()->SetId (j,id);
	  
	     //insert field
	     STATE* sl = (STATE*)left_state(p);
	     impulses->InsertTuple(id,sl->collsnImpulse);
	     if (sl->collsn_num > 0)
	         has_collsn = true;
	 }
         cells->InsertNextCell(cell);
	 if (has_collsn) colors->InsertNextTupleValue(red);
	 else colors->InsertNextTupleValue(green);
  }
  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();

  polydata->SetPolys(cells);
  polydata->GetCellData()->SetScalars(colors);
  polydata->SetPoints(pts);
  polydata->GetPointData()->SetVectors(impulses);

  std::cout<<"vtk plot: "<<polydata->GetNumberOfPoints()<<" points"<<std::endl;
  std::cout<<"vtk plot: "<<polydata->GetNumberOfCells()<<" cells"<<std::endl;
  // Write the file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(fname);
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(polydata);
#else
  writer->SetInputData(polydata);
#endif
  // Optional - set the mode. The default is binary.
  writer->SetDataModeToBinary();
  //writer->SetDataModeToAscii();
  writer->Write();
  
  return;
}
