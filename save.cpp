#include "include/traj.h"

void save_files(int i_time, unsigned int n_space_div[3], float posL[3], float dd[3], double t,
                float np[2][n_space_divz][n_space_divy][n_space_divx], float currentj[2][3][n_space_divz][n_space_divy][n_space_divx],
                float V[n_space_divz][n_space_divy][n_space_divx],
                float E[3][n_space_divz][n_space_divy][n_space_divx], float B[3][n_space_divz][n_space_divy][n_space_divx],
                float KE[2][n_output_part], float posp[2][n_output_part][3])
{
#ifdef printDensity
  save_vti_c("Ne", i_time, n_space_div, posL, dd, n_cells, 1, t, &np[0], "Float32", sizeof(float));
  save_vti_c("je", i_time, n_space_div, posL, dd, n_cells, 3, t, currentj[0], "Float32", sizeof(float));
#endif
#ifdef printV
  save_vti_c("V", i_time, n_space_div, posL, dd, n_cells, 1, t, V, "Float32", sizeof(float));
#endif
#ifdef printE
  save_vti_c("E", i_time, n_space_div, posL, dd, n_cells, 3, t, E, "Float32", sizeof(float));
#endif
#ifdef printB
  save_vti_c("B", i_time, n_space_div, posL, dd, n_cells, 3, t, B, "Float32", sizeof(float));
#endif
#ifdef printParticles
  //  save_vtp("e", i_time, n_output_part, 1, t, (reinterpret_cast<const char *>(&KE[0][0])), (reinterpret_cast<const char *>(&posp[0][0][0])));
  // save_vtp("d", i_time, n_output_part, 1, t, (reinterpret_cast<const char *>(&KE[1][0])), (reinterpret_cast<const char *>(&posp[1][0][0])));
  save_vtp("e", i_time, n_output_part, 0, t, KE, posp);
  save_vtp("d", i_time, n_output_part, 1, t, KE, posp);
#endif
}
void save_hist(int i, double t, int npart, int mp[2], float dt[2], float pos0x[2][n_partd], float pos0y[2][n_partd], float pos0z[2][n_partd], float pos1x[2][n_partd], float pos1y[2][n_partd], float pos1z[2][n_partd])
{
  // Create the vtkTable object
  vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
  table->SetNumberOfRows(Hist_n);
  // Create the time array
  // vtkSmartPointer<vtkDoubleArray> timeArray = vtkSmartPointer<vtkDoubleArray>::New();
  // timeArray->SetName("Time");
  // timeArray->InsertNextValue(t);
  // table->AddColumn(timeArray);

  // Create the histogram arrays
  vtkSmartPointer<vtkDoubleArray> energyArray = vtkSmartPointer<vtkDoubleArray>::New();
  energyArray->SetName("Energy(eV)");
  vtkSmartPointer<vtkDoubleArray> electronHistArray = vtkSmartPointer<vtkDoubleArray>::New();
  electronHistArray->SetName("Electron KE Histogram");
  vtkSmartPointer<vtkDoubleArray> ionHistArray = vtkSmartPointer<vtkDoubleArray>::New();
  ionHistArray->SetName("Ion KE Histogram");

  double KEhist[2][Hist_n];
  memset(KEhist, 0, sizeof(KEhist));
  for (int p = 0; p < 2; ++p)
    for (int i = 0; i < npart; ++i)
    {
      float dx = pos1x[p][i] - pos0x[p][i];
      float dy = pos1z[p][i] - pos0y[p][i];
      float dz = pos1z[p][i] - pos0z[p][i];
      int index = (int)trunc((0.5 * (float)mp[p] * (dx * dx + dy * dy + dz * dz) / (e_charge_mass * dt[p] * dt[p] * (float)Hist_n)) /(float) Hist_max);
      if (index > Hist_n)
        index = Hist_n - 1;
      if (index < 0)
        cout << "error index<0" << endl;
      KEhist[p][index]++;
    }

  // Add the histogram values to the arrays
  for (int i = 0; i < Hist_n; ++i)
  {
    energyArray->InsertNextValue((double)(i * Hist_max) / (double)(Hist_n));
    electronHistArray->InsertNextValue(KEhist[0][i] + 1);
    ionHistArray->InsertNextValue(KEhist[1][i] + 1);
  }

  // Add the histogram arrays to the table
  table->AddColumn(energyArray);
  table->AddColumn(electronHistArray);
  table->AddColumn(ionHistArray);

  // Write the table to a file
  vtkSmartPointer<vtkDelimitedTextWriter> writer = vtkSmartPointer<vtkDelimitedTextWriter>::New();
  writer->SetFileName((outpath + "KEhist_" + to_string(i) + ".csv").c_str());
  writer->SetInputData(table);
  writer->Write();
}

/*
void save_hist(double t, int npart, int mp[2], float dt[2], float pos0x[2][n_partd], float pos0y[2][n_partd], float pos0z[2][n_partd], float pos1x[2][n_partd], float pos1y[2][n_partd], float pos1z[2][n_partd])
{
  ofstream Histe_file, Histd_file;
  float KEhist[2][Hist_n];
  memset(KEhist, 0, sizeof(KEhist));
  static int first = 1;
  if (first)
  {
    Histe_file.open(outpath + "Histe.csv");
    Histd_file.open(outpath + "Histd.csv");
    Histe_file << "t";
    Histd_file << "t";
    for (int i = 0; i < Hist_n; ++i)
    {
      // Histe_file << ",Ee" << i;
      // Histd_file << ",Ed" << i;
      Histe_file << "," << i;
      Histd_file << "," << i;
    }
    Histe_file << endl;
    Histd_file << endl;
    Histe_file.close();
    Histd_file.close();
    first = 0;
  }
  Histe_file.open(outpath + "Histe.csv", std::ios_base::app);
  Histd_file.open(outpath + "Histd.csv", std::ios_base::app);
  Histe_file << t;
  Histd_file << t;

  for (int p = 0; p < 2; ++p)
    for (int i = 0; i < npart; ++i)
    {
      float dx = pos1x[p][i] - pos0x[p][i];
      float dy = pos1z[p][i] - pos0y[p][i];
      float dz = pos1z[p][i] - pos0z[p][i];
      int index = (int)trunc((0.5 * (float)mp[p] * (dx * dx + dy * dy + dz * dz) / (e_charge_mass * dt[p] * dt[p] * (float)Hist_n)) / Hist_max);
      if (index > Hist_n)
        index = Hist_n - 1;
      if (index < 0)
        cout << "error index<0" << endl;
      KEhist[p][index]++;
    }
  for (int i = 0; i < Hist_n; ++i)
  {
    Histe_file << "," << KEhist[0][i];
    Histd_file << "," << KEhist[1][i];
  }
  Histe_file << endl;
  Histd_file << endl;
  Histe_file.close();
  Histd_file.close();
}
*/
/**
 * This corrects the order of dimensions for view in paraview, as opposed to save_vti which prints the raw data.
 */
void save_vti_c(string filename, int i,
                unsigned int n_space_div[3], float posl[3], float dd[3], uint64_t num, int ncomponents, double t,
                float data1[][n_space_divz][n_space_divy][n_space_divz], string typeofdata, int bytesperdata)
{
  if (ncomponents > 3)
  {
    cout << "Error: Cannot write file " << filename << " - too many components" << endl;
    return;
  }

  int xi = (n_space_divx - 1) / maxcells + 1;
  int yj = (n_space_divy - 1) / maxcells + 1;
  int zk = (n_space_divz - 1) / maxcells + 1;
  int nx = n_space_div[0] / xi;
  int ny = n_space_div[1] / yj;
  int nz = n_space_div[2] / zk;
  // vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New(); // Create the vtkImageData object
  imageData->SetDimensions(nx, ny, nz);                                           // Set the dimensions of the image data
  imageData->SetSpacing(dd[0] * xi, dd[1] * yj, dd[2] * zk);
  imageData->SetOrigin(posl[0], posl[1], posl[2]); // Set the origin of the image data
  imageData->AllocateScalars(VTK_FLOAT, ncomponents);
  imageData->GetPointData()->GetScalars()->SetName(filename.c_str());
  float *data2 = static_cast<float *>(imageData->GetScalarPointer()); // Get a pointer to the density field array

  for (int k = 0; k < nz; ++k)
  {
    for (int j = 0; j < ny; ++j)
    {
      for (int i = 0; i < nx; ++i)
      {
        for (int c = 0; c < ncomponents; ++c)
          data2[(k * ny + j) * nx * ncomponents + i * ncomponents + c] = data1[c][k * zk][j * yj][i * xi];
      }
    }
  }

  vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New(); // Create the vtkXMLImageDataWriter object
  writer->SetFileName((outpath + filename + "_" + to_string(i) + ".vti").c_str());               // Set the output file name                                                                     // Set the time value
  writer->SetDataModeToBinary();
  // writer->SetCompressorTypeToLZ4();
  writer->SetCompressorTypeToZLib(); // Enable compression
  writer->SetCompressionLevel(9);    // Set the level of compression (0-9)
  writer->SetInputData(imageData);   // Set the input image data
  writer->Write();                   // Write the output file
}

void save_vtp(string filename, int i, uint64_t num, int n, double t, float data[2][n_output_part], float points1[2][n_output_part][3])
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkFloatArray> kineticEnergy = vtkSmartPointer<vtkFloatArray>::New();
  kineticEnergy->SetName("KE");
  for (int i = 0; i < num; ++i)
  {
    points->InsertNextPoint(points1[n][i][0], points1[n][i][1], points1[n][i][2]);
    kineticEnergy->InsertNextValue(data[n][i]);
  }
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New(); // Create the VTK poly data object
  polyData->SetPoints(points);
  polyData->GetPointData()->AddArray(kineticEnergy);
  // Write the output file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName((outpath + filename + "_" + to_string(i) + ".vtp").c_str());
  writer->SetDataModeToBinary();
  writer->SetCompressorTypeToZLib(); // Enable compression
  writer->SetCompressionLevel(9);    // Set the level of compression (0-9)
  writer->SetInputData(polyData);
  writer->Write();
}

void save_vti_c2(string filename, int i,
                 unsigned int n_space_div[3], float posl[3], float dd[3], uint64_t num, int ncomponents, double t,
                 float data1[3][n_space_divz2][n_space_divy2][n_space_divz2], string typeofdata, int bytesperdata)
{ // TODO: fix dimension sizes
  auto *data = new float[n_space_divz2][n_space_divy2][n_space_divx2][3];
  for (int k = 0; k < n_space_divz2; ++k)
  {
    for (int j = 0; j < n_space_divy2; ++j)
    {
      for (int i = 0; i < n_space_divx2; ++i)
      {
        for (int c = 0; c < 3; ++c)
          data[k][j][i][c] = data1[c][k][j][i];
      }
    }
  }

  //  cout <<std::filesystem::temp_directory_path().string()+"/out/"+filename+"_"+to_string(i)+".vti";
  std::ofstream os(outpath + filename + "_" + to_string(i) + ".vti", std::ios::binary | std::ios::out);
  os << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\"> \n ";
  os << "<ImageData WholeExtent=\"0 ";
  os << to_string(n_space_div[0] - 1) + " 0 " + to_string(n_space_div[1] - 1) + " 0 " + to_string(n_space_div[2] - 1) + "\" ";
  os << "Origin=\"" + to_string(posl[0]) + " " + to_string(posl[1]) + " " + to_string(posl[2]) + "\"";
  os << " Spacing=\"" + to_string(dd[0]) + " " + to_string(dd[1]) + " " + to_string(dd[2]) + "\" ";
  os << "Direction=\"1 0 0 0 1 0 0 0 1\"> \n";
  os << "<FieldData>\n";
  os << "<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" offset=\"0\"/>\n";
  os << "</FieldData>\n";
  os << "<Piece Extent=\"0 ";
  os << to_string(n_space_div[0] - 1) + " 0 " + to_string(n_space_div[1] - 1) + " 0 " + to_string(n_space_div[2] - 1) + "\">\n";
  os << "<PointData Scalars=\"" + filename + "\">\n";
  os << "<DataArray type=\"" + typeofdata + "\" Name=\"" + filename + "\" NumberOfComponents=\"" + to_string(ncomponents) + "\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" offset=\"16\" />\n";
  os << "  </PointData>\n";
  os << "<CellData>\n";
  os << "  </CellData>\n";
  os << "</Piece>\n";
  os << "</ImageData>\n";
  os << "<AppendedData encoding=\"raw\">_";
  uint64_t num1 = 8;
  os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
  // single time double
  os.write(reinterpret_cast<const char *>(&t), std::streamsize(sizeof(double)));
  num1 = num * ncomponents * bytesperdata;
  os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
  // data
  os.write((reinterpret_cast<const char *>(data)), std::streamsize(num * ncomponents * bytesperdata));
  os << "</AppendedData>\n";
  os << "</VTKFile>";
  os.close();
  delete[] data;
}

void save_pvd(string filename, int ndatapoints)
{
  std::ofstream os(outpath + filename + ".pvd", std::ios::binary | std::ios::out);
  os << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  os << "<Collection>\n";
  for (int i = 0; i < ndatapoints; i++)
  {
    os << "<DataSet timestep=\"" + to_string(i) + "\" part=\"0\" file=\"" + filename + "_" + to_string(i) + ".vti\"/>\n";
  }
  os << " </Collection>\n";
  os << "</VTKFile>\n";
  os.close();
}
/*
void save_vti(string filename, int i,
              unsigned int n_space_div[3], float posl[3], float dd[3], uint64_t num, int ncomponents, double t,
              float data[n_space_divz][n_space_divy][n_space_divx], string typeofdata, int bytesperdata)
{
  int nx = n_space_div[0];
  int ny = n_space_div[1];
  int nz = n_space_div[2];
  vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New(); // Create the vtkImageData object
  imageData->SetDimensions(nx, ny, nz);                                           // Set the dimensions of the image data
  imageData->SetSpacing(dd[0], dd[1], dd[2]);
  imageData->SetOrigin(posl[0], posl[1], posl[2]); // Set the origin of the image data
  imageData->AllocateScalars(VTK_FLOAT, 1);
  float *densityArray = static_cast<float *>(imageData->GetScalarPointer()); // Get a pointer to the density field array
                                                                             //  imageData->GetFieldData()->SetObjectName("aaa");
  imageData->GetPointData()->GetScalars()->SetName(filename.c_str());
  // Fill the density field array with your data
  for (int x = 0; x < nx; x++)
  {
    for (int y = 0; y < ny; y++)
    {
      for (int z = 0; z < nz; z++)
      {
        densityArray[(z * ny + y) * nx + x] = data[x][y][z]; // Set the density value at this point
      }
    }
  }

  vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New(); // Create the vtkXMLImageDataWriter object
  writer->SetFileName((outpath + filename + "_" + to_string(i) + ".vti").c_str());               // Set the output file name
  writer->SetCompressorTypeToZLib();                                                             // Enable compression
  writer->SetCompressionLevel(9);                                                                // Set the level of compression (0-9)
  writer->SetInputData(imageData);                                                               // Set the input image data
  writer->Write();                                                                               // Write the output file
}
*/
/*
//  cout <<std::filesystem::temp_directory_path().string()+"/out/"+filename+"_"+to_string(i)+".vti";
std::ofstream os(outpath + filename + "_" + to_string(i) + ".vti", std::ios::binary | std::ios::out);
os << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\"> \n ";
os << "<ImageData WholeExtent=\"0 ";
os << to_string(n_space_div[0] / xi - 1) + " 0 " + to_string(n_space_div[1] / yj - 1) + " 0 " + to_string(n_space_div[2] / zk - 1) + "\" ";
os << "Origin=\"" + to_string(posl[0]) + " " + to_string(posl[1]) + " " + to_string(posl[2]) + "\"";
os << " Spacing=\"" + to_string(dd[0] * xi) + " " + to_string(dd[1] * yj) + " " + to_string(dd[2] * zk) + "\" ";
os << "Direction=\"1 0 0 0 1 0 0 0 1\"> \n";
os << "<FieldData>\n";
os << "<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" offset=\"0\"/>\n";
os << "</FieldData>\n";
os << "<Piece Extent=\"0 ";
os << to_string(n_space_div[0] / xi - 1) + " 0 " + to_string(n_space_div[1] / yj - 1) + " 0 " + to_string(n_space_div[2] / zk - 1) + "\">\n";
os << "<PointData Scalars=\"" + filename + "\">\n";
os << "<DataArray type=\"" + typeofdata + "\" Name=\"" + filename + "\" NumberOfComponents=\"" + to_string(ncomponents) + "\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" offset=\"16\" />\n";
os << "  </PointData>\n";
os << "<CellData>\n";
os << "  </CellData>\n";
os << "</Piece>\n";
os << "</ImageData>\n";
os << "<AppendedData encoding=\"raw\">_";
uint64_t num1 = 8;
os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
// single time double
os.write(reinterpret_cast<const char *>(&t), std::streamsize(sizeof(double)));
num1 = num * ncomponents * bytesperdata / xi / yj / zk;
os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
// data
os.write((reinterpret_cast<const char *>(data)), std::streamsize(num * ncomponents * bytesperdata / xi / yj / zk));
os << "</AppendedData>\n";
os << "</VTKFile>";
os.close();
delete[] data;
*/
/*
std::ofstream os(outpath + filename + "_" + to_string(i) + ".vtp", std::ios::binary | std::ios::out);
os << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
os << "<PolyData>\n <FieldData>\n  <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" offset=\"0\"/>\n</FieldData>\n";
os << "<Piece NumberOfPoints=\"" + to_string(num);
os << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\" >\n";
os << " <PointData>\n";
os << "  <DataArray type=\"Float32\" Name=\" KE\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" ";
os << "offset=\"16\"/>\n";
os << " </PointData>\n  <CellData>\n </CellData>\n";
os << "<Points>";
os << "    <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" RangeMin=\"0.0\" RangeMax=\"0.0\" offset=\"" + to_string(num * sizeof(float) * ncomponents + 24) + "\"/>\n";
os << "  </Points>\n";
os << "</Piece>\n";
os << "</PolyData>\n";
os << "<AppendedData encoding=\"raw\">\n_";
uint64_t num1 = 8;
os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
os.write(reinterpret_cast<const char *>(&t), std::streamsize(sizeof(double)));
num1 = num * ncomponents * sizeof(float);
os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
os.write(data, std::streamsize(num * ncomponents * sizeof(float)));
num1 = num * sizeof(float) * 3;
os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
os.write(points, std::streamsize(num * 3 * sizeof(float)));
os << "\n</AppendedData>\n";
os << "</VTKFile>";
os.close();
*/