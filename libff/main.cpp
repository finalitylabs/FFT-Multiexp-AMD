#include <cstdio>
#include <vector>

#include <libff/algebra/curves/mnt753/mnt4753/mnt4753_pp.hpp>
#include <libff/algebra/curves/mnt753/mnt6753/mnt6753_pp.hpp>

// #define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.h>
#define DATA_SIZE (131072)
#define limbs_per_elem (12)
#include <chrono> 

using namespace std::chrono; 
using namespace std;

#include <typeinfo>
#include <unistd.h>
#include <string.h>
#include <errno.h>

char *getcwd(char *buf, size_t size);

using namespace libff;

// mnt4 montgomery
Fq<mnt4753_pp> read_mnt4_fq(FILE* input) {
  Fq<mnt4753_pp> x;
  fread((void *) x.mont_repr.data, libff::mnt4753_q_limbs * sizeof(mp_size_t), 1, input);
  return x;
}

void write_mnt4_fq(FILE* output, Fq<mnt4753_pp> x) {
  fwrite((void *) x.mont_repr.data, libff::mnt4753_q_limbs * sizeof(mp_size_t), 1, output);
}

// mnt6 montgomery
Fq<mnt6753_pp> read_mnt6_fq(FILE* input) {
  Fq<mnt6753_pp> x;
  fread((void *) x.mont_repr.data, libff::mnt6753_q_limbs * sizeof(mp_size_t), 1, input);
  return x;
}

void write_mnt6_fq(FILE* output, Fq<mnt6753_pp> x) {
  fwrite((void *) x.mont_repr.data, libff::mnt6753_q_limbs * sizeof(mp_size_t), 1, output);
}

// mnt4 fq2 montgomery
Fqe<mnt4753_pp> read_mnt4_fq2(FILE* input) {
  Fq<mnt4753_pp> c0 = read_mnt4_fq(input);
  Fq<mnt4753_pp> c1 = read_mnt4_fq(input);
  return Fqe<mnt4753_pp>(c0, c1);
}

void write_mnt4_fq2(FILE* output, Fqe<mnt4753_pp> x) {
  write_mnt4_fq(output, x.c0);
  write_mnt4_fq(output, x.c1);
}

// mnt6 fq3 montgomery
Fqe<mnt6753_pp> read_mnt6_fq3(FILE* input) {
  Fq<mnt6753_pp> c0 = read_mnt6_fq(input);
  Fq<mnt6753_pp> c1 = read_mnt6_fq(input);
  Fq<mnt6753_pp> c2 = read_mnt6_fq(input);
  return Fqe<mnt6753_pp>(c0, c1, c2);
}

void write_mnt6_fq3(FILE* output, Fqe<mnt6753_pp> x) {
  write_mnt6_fq(output, x.c0);
  write_mnt6_fq(output, x.c1);
  write_mnt6_fq(output, x.c2);
}

// mnt4 groups montgomery
G1<mnt4753_pp> read_mnt4_g1(FILE* input) {
  Fq<mnt4753_pp> x = read_mnt4_fq(input);
  Fq<mnt4753_pp> y = read_mnt4_fq(input);
  return G1<mnt4753_pp>(x, y, Fq<mnt4753_pp>::one());
}

void write_mnt4_g1(FILE* output, G1<mnt4753_pp> g) {
  g.to_affine_coordinates();
  write_mnt4_fq(output, g.X());
  write_mnt4_fq(output, g.Y());
}

G2<mnt4753_pp> read_mnt4_g2(FILE* input) {
  Fqe<mnt4753_pp> x = read_mnt4_fq2(input);
  Fqe<mnt4753_pp> y = read_mnt4_fq2(input);
  return G2<mnt4753_pp>(x, y, Fqe<mnt4753_pp>::one());
}

void write_mnt4_g2(FILE* output, G2<mnt4753_pp> g) {
  g.to_affine_coordinates();
  write_mnt4_fq2(output, g.X());
  write_mnt4_fq2(output, g.Y());
}

// mnt6 groups montgomery
G1<mnt6753_pp> read_mnt6_g1(FILE* input) {
  Fq<mnt6753_pp> x = read_mnt6_fq(input);
  Fq<mnt6753_pp> y = read_mnt6_fq(input);
  return G1<mnt6753_pp>(x, y, Fq<mnt6753_pp>::one());
}

void write_mnt6_g1(FILE* output, G1<mnt6753_pp> g) {
  g.to_affine_coordinates();
  write_mnt6_fq(output, g.X());
  write_mnt6_fq(output, g.Y());
}

G2<mnt6753_pp> read_mnt6_g2(FILE* input) {
  Fqe<mnt6753_pp> x = read_mnt6_fq3(input);
  Fqe<mnt6753_pp> y = read_mnt6_fq3(input);
  return G2<mnt6753_pp>(x, y, Fqe<mnt6753_pp>::one());
}

void write_mnt6_g2(FILE* output, G2<mnt6753_pp> g) {
  g.to_affine_coordinates();
  write_mnt6_fq3(output, g.X());
  write_mnt6_fq3(output, g.Y());
}

// mnt4 ordinary
Fq<mnt4753_pp> read_mnt4_fq_numeral(FILE* input) {
  // bigint<mnt4753_q_limbs> n;
  Fq<mnt4753_pp> x;
  fread((void *) x.mont_repr.data, libff::mnt4753_q_limbs * sizeof(mp_size_t), 1, input);
  auto b = Fq<mnt4753_pp>(x.mont_repr);
  return b;
}

void write_mnt4_fq_numeral(FILE* output, Fq<mnt4753_pp> x) {
  auto out_numeral = x.as_bigint();
  fwrite((void *) out_numeral.data, libff::mnt4753_q_limbs * sizeof(mp_size_t), 1, output);
}

// mnt6 ordinary
Fq<mnt6753_pp> read_mnt6_fq_numeral(FILE* input) {
  // bigint<mnt4753_q_limbs> n;
  Fq<mnt6753_pp> x;
  fread((void *) x.mont_repr.data, libff::mnt6753_q_limbs * sizeof(mp_size_t), 1, input);
  auto b = Fq<mnt6753_pp>(x.mont_repr);
  return b;
}

void write_mnt6_fq_numeral(FILE* output, Fq<mnt6753_pp> x) {
  auto out_numeral = x.as_bigint();
  fwrite((void *) out_numeral.data, libff::mnt6753_q_limbs * sizeof(mp_size_t), 1, output);
}

// mnt4 fq2 ordinary
Fqe<mnt4753_pp> read_mnt4_fq2_numeral(FILE* input) {
  Fq<mnt4753_pp> c0 = read_mnt4_fq_numeral(input);
  Fq<mnt4753_pp> c1 = read_mnt4_fq_numeral(input);
  return Fqe<mnt4753_pp>(c0, c1);
}

void write_mnt4_fq2_numeral(FILE* output, Fqe<mnt4753_pp> x) {
  write_mnt4_fq_numeral(output, x.c0);
  write_mnt4_fq_numeral(output, x.c1);
}

// mnt6 fq3 ordinary
Fqe<mnt6753_pp> read_mnt6_fq3_numeral(FILE* input) {
  Fq<mnt6753_pp> c0 = read_mnt6_fq_numeral(input);
  Fq<mnt6753_pp> c1 = read_mnt6_fq_numeral(input);
  Fq<mnt6753_pp> c2 = read_mnt6_fq_numeral(input);
  return Fqe<mnt6753_pp>(c0, c1, c2);
}

void write_mnt6_fq3_numeral(FILE* output, Fqe<mnt6753_pp> x) {
  write_mnt6_fq_numeral(output, x.c0);
  write_mnt6_fq_numeral(output, x.c1);
  write_mnt6_fq_numeral(output, x.c2);
}

// mnt4 groups ordinary
G1<mnt4753_pp> read_mnt4_g1_numeral(FILE* input) {
  Fq<mnt4753_pp> x = read_mnt4_fq_numeral(input);
  Fq<mnt4753_pp> y = read_mnt4_fq_numeral(input);
  return G1<mnt4753_pp>(x, y, Fq<mnt4753_pp>::one());
}

void write_mnt4_g1_numeral(FILE* output, G1<mnt4753_pp> g) {
  g.to_affine_coordinates();
  write_mnt4_fq_numeral(output, g.X());
  write_mnt4_fq_numeral(output, g.Y());
}

G2<mnt4753_pp> read_mnt4_g2_numeral(FILE* input) {
  Fqe<mnt4753_pp> x = read_mnt4_fq2_numeral(input);
  Fqe<mnt4753_pp> y = read_mnt4_fq2_numeral(input);
  return G2<mnt4753_pp>(x, y, Fqe<mnt4753_pp>::one());
}

void write_mnt4_g2_numeral(FILE* output, G2<mnt4753_pp> g) {
  g.to_affine_coordinates();
  write_mnt4_fq2_numeral(output, g.X());
  write_mnt4_fq2_numeral(output, g.Y());
}

// mnt6 groups ordinary
G1<mnt6753_pp> read_mnt6_g1_numeral(FILE* input) {
  Fq<mnt6753_pp> x = read_mnt6_fq_numeral(input);
  Fq<mnt6753_pp> y = read_mnt6_fq_numeral(input);
  return G1<mnt6753_pp>(x, y, Fq<mnt6753_pp>::one());
}

void write_mnt6_g1_numeral(FILE* output, G1<mnt6753_pp> g) {
  g.to_affine_coordinates();
  write_mnt6_fq_numeral(output, g.X());
  write_mnt6_fq_numeral(output, g.Y());
}

G2<mnt6753_pp> read_mnt6_g2_numeral(FILE* input) {
  Fqe<mnt6753_pp> x = read_mnt6_fq3_numeral(input);
  Fqe<mnt6753_pp> y = read_mnt6_fq3_numeral(input);
  return G2<mnt6753_pp>(x, y, Fqe<mnt6753_pp>::one());
}

void write_mnt6_g2_numeral(FILE* output, G2<mnt6753_pp> g) {
  g.to_affine_coordinates();
  write_mnt6_fq3_numeral(output, g.X());
  write_mnt6_fq3_numeral(output, g.Y());
}

// The actual code for doing the group operations lives in
// libff/algebra/curves/mnt753/mnt4753/mnt4753_g1.tcc
// libff/algebra/curves/mnt753/mnt4753/mnt4753_g2.tcc
// libff/algebra/curves/mnt753/mnt6753/mnt6753_g1.tcc
// libff/algebra/curves/mnt753/mnt6753/mnt6753_g2.tcc
int main(int argc, char *argv[])
{
    // argv should be
    // { "main", "compute" or "compute-numeral", inputs, outputs }

    mnt4753_pp::init_public_params();
    mnt6753_pp::init_public_params();

    size_t n;

    auto is_numeral = strcmp(argv[1], "compute-numeral") == 0;
    auto inputs = fopen(argv[2], "r");
    auto outputs = fopen(argv[3], "w");

    // mnt4
    auto read_mnt4 = read_mnt4_fq;
    auto write_mnt4 = write_mnt4_fq;
    auto read_mnt4_q2 = read_mnt4_fq2;
    auto write_mnt4_q2 = write_mnt4_fq2;
    // mnt6
    auto read_mnt6 = read_mnt6_fq;
    auto write_mnt6 = write_mnt6_fq;
    auto read_mnt6_q3 = read_mnt6_fq3;
    auto write_mnt6_q3 = write_mnt6_fq3;
    // mnt4 groups
    auto write_mnt4g1 = write_mnt4_g1;
    auto read_mnt4g1 = read_mnt4_g1;
    auto write_mnt4g2 = write_mnt4_g2;
    auto read_mnt4g2 = read_mnt4_g2;
    // mnt6 groups
    auto read_mnt6g1 = read_mnt6_g1;
    auto write_mnt6g1 = write_mnt6_g1;
    auto read_mnt6g2 = read_mnt6_g2;
    auto write_mnt6g2 = write_mnt6_g2;

    if (is_numeral) {
      // mnt4
      read_mnt4 = read_mnt4_fq_numeral;
      write_mnt4 = write_mnt4_fq_numeral;
      read_mnt4_q2 = read_mnt4_fq2_numeral;
      write_mnt4_q2 = write_mnt4_fq2_numeral;
      // mnt6
      read_mnt6 = read_mnt6_fq_numeral;
      write_mnt6 = write_mnt6_fq_numeral;
      read_mnt6_q3 = read_mnt6_fq3_numeral;
      write_mnt6_q3 = write_mnt6_fq3_numeral;
      // mnt4 groups
      write_mnt4g1 = write_mnt4_g1_numeral;
      read_mnt4g1 = read_mnt4_g1_numeral;
      write_mnt4g2 = write_mnt4_g2_numeral;
      read_mnt4g2 = read_mnt4_g2_numeral;
      // mnt6 groups
      read_mnt6g1 = read_mnt6_g1_numeral;
      write_mnt6g1 = write_mnt6_g1_numeral;
      read_mnt6g2 = read_mnt6_g2_numeral;
      write_mnt6g2 = write_mnt6_g2_numeral;
    }


    while (true) {
      size_t elts_read = fread((void *) &n, sizeof(size_t), 1, inputs);
      if (elts_read == 0) { break; }

      // Read input
      std::vector<G1<mnt4753_pp>> g4_1;
      for (size_t i = 0; i < n; ++i) { g4_1.emplace_back(read_mnt4g1(inputs)); }

      std::vector<G2<mnt4753_pp>> g4_2;
      for (size_t i = 0; i < n; ++i) { g4_2.emplace_back(read_mnt4g2(inputs)); }

      std::vector<G1<mnt6753_pp>> g6_1;
      for (size_t i = 0; i < n; ++i) { g6_1.emplace_back(read_mnt6g1(inputs)); }

      std::vector<G2<mnt6753_pp>> g6_2;
      for (size_t i = 0; i < n; ++i) { g6_2.emplace_back(read_mnt6g2(inputs)); }

      // Perform the computation
      G1<mnt4753_pp> h4_1 = G1<mnt4753_pp>::zero();
      printf("g1_h4 zero x\n");
      h4_1.print();

      for (size_t i = 0; i < n; ++i) { h4_1 = h4_1 + g4_1[i]; }

      G2<mnt4753_pp> h4_2 = G2<mnt4753_pp>::zero();
      for (size_t i = 0; i < n; ++i) { h4_2 = h4_2 + g4_2[i]; }

      G1<mnt6753_pp> h6_1 = G1<mnt6753_pp>::zero();
      for (size_t i = 0; i < n; ++i) { h6_1 = h6_1 + g6_1[i]; }

      G2<mnt6753_pp> h6_2 = G2<mnt6753_pp>::zero();
      for (size_t i = 0; i < n; ++i) { h6_2 = h6_2 + g6_2[i]; }

      // Write output
      write_mnt4g1(outputs, h4_1);
      write_mnt4g2(outputs, h4_2);
      write_mnt6g1(outputs, h6_1);
      write_mnt6g2(outputs, h6_2);

      printf("g1_h4  x\n");
      h4_1.print();

      printf("g1_g4 y\n");
      g4_1[498].print();

      // OPENCL START

      char cwd[1024];
      if (getcwd(cwd, sizeof(cwd)) != NULL) {
         // printf("Current working dir: %s\n", cwd);
      } else {
         perror("getcwd() error");
         return 1;
      }

      FILE *fp;
      char *source_str;
      size_t source_size, program_size;

      fp = fopen("kernels/ec.cl", "r");
      if (!fp) {
          fprintf(stderr, "could not open program file\n");
          exit(1);
      }

      char* program_source_code;
      size_t program_source_code_size;
      program_source_code = (char*)malloc(400000);
      program_source_code_size = fread(program_source_code, 1, 400000, fp);
      fclose(fp);

      int err;                            // error code returned from api calls
      char name[128];
        
      G1<mnt4753_pp>* data_x = new G1<mnt4753_pp>[1];              // original data set given to device
      G1<mnt4753_pp>* data_y = new G1<mnt4753_pp>[n];              // original data set given to device
      G1<mnt4753_pp>* results = new G1<mnt4753_pp>[1];           // results returned from device
      unsigned int correct;               // number of correct results returned

      size_t global;                      // global domain size for our calculation
      size_t local;                       // local domain size for our calculation
      cl_device_id device_id;             // compute device id 
      cl_context context;                 // compute context
      cl_command_queue commands;          // compute command queue
      cl_program program;                 // compute program
      cl_kernel kernel;                   // compute kernel
      cl_event event;                     // timing
      cl_ulong time_start;
      cl_ulong time_end;

      
      cl_mem input_x;                       // device memory used for the input array
      cl_mem input_y;                       // device memory used for the input array
      cl_mem output;                       // device memory used for the input array
      // Fill our data set with field inputs from param gen
      //
      unsigned int count = n;
      mp_size_t num = 1;

      //memcpy(&data_x[0], &h4_1, sizeof(G1<mnt4753_pp>));
      data_x[0] = G1<mnt4753_pp>::zero();
      printf("count %u\n", n);
      data_x[0].print_coordinates();

      for(int i = 0; i < count; i++) {
        memcpy(&data_y[i], &g4_1[i], sizeof(G1<mnt4753_pp>));
      }
      
      // Connect to a compute device
      //

      /* get platform number of OpenCL */
      cl_uint  num_platforms = 0;
      clGetPlatformIDs (0, NULL, &num_platforms);
      printf("num_platforms: %d\n", (int)num_platforms);

      /* allocate a segment of mem space, so as to store supported platform info */
      cl_platform_id *platforms = (cl_platform_id *) malloc (num_platforms * sizeof (cl_platform_id));

      /* get platform info */
      clGetPlatformIDs (num_platforms, platforms, NULL);

      /* get device number on platform */
      cl_uint num_devices = 0;
      clGetDeviceIDs (platforms[0], CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices);
      printf("num_devices: %d\n", (int)num_devices);

      /* allocate a segment of mem space, to store device info, supported by platform */
      cl_device_id *devices;
      devices = (cl_device_id *) malloc (num_devices * sizeof (cl_device_id));

      /* get device info */
      clGetDeviceIDs (platforms[0], CL_DEVICE_TYPE_GPU, num_devices, devices, NULL);

      // int gpu = 1;
      // err = clGetDeviceIDs(NULL, gpu ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 1, &device_id, NULL);
      // if (err != CL_SUCCESS)
      // {
      //     printf("Error: Failed to create a device group!\n");
      //     return EXIT_FAILURE;
      // }

      printf("Device id: %u\n", devices[0]);

      clGetDeviceInfo(devices[0], CL_DEVICE_NAME, 128, name, NULL);
      fprintf(stdout, "Created a dispatch queue using the %s\n", name);

      // Create a compute context 
      //
      printf("creating context\n");
      context = clCreateContext(0, num_devices, devices, NULL, NULL, &err);
      if (!context)
      {
          printf("Error: Failed to create a compute context!\n");
          return EXIT_FAILURE;
      }

      // Create a command commands
      //
      commands = clCreateCommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE, &err);
      if (!commands)
      {
          printf("Error: Failed to create a command commands!\n");
          return EXIT_FAILURE;
      }

      // Create the compute program from the source buffer
      //
      program = clCreateProgramWithSource(context, 1, (const char **) &program_source_code, &program_source_code_size, &err);
      if (!program)
      {
          printf("Error: Failed to create compute program!\n");
          return EXIT_FAILURE;
      }

      // Build the program executable
      //
      printf("building program\n");
      char options[] = "-cl-opt-disable";
      err = clBuildProgram(program, num_devices, devices, NULL, NULL, NULL);
      if (err != CL_SUCCESS)
      {
          size_t len;
          char buffer[2048];
          //std::cerr << getErrorString(err) << std::endl;
          printf("Error: Failed to build program executable!\n");
          printf ("Message: %s\n",strerror(err));
          clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
          exit(1);
      }

      // Create the compute kernel in the program we wish to run
      //
      kernel = clCreateKernel(program, "add_G1", &err);
      if (!kernel || err != CL_SUCCESS)
      {
          printf("Error: Failed to create compute kernel!\n");
          exit(1);
      }

      // Create the input and output arrays in device memory for our calculation
      //
      printf("creating buffer\n");
      input_x = clCreateBuffer(context,  CL_MEM_READ_ONLY,  sizeof(G1<mnt4753_pp>), NULL, NULL);
      input_y = clCreateBuffer(context,  CL_MEM_READ_ONLY,  sizeof(G1<mnt4753_pp>) * count, NULL, NULL);
      output = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(G1<mnt4753_pp>), NULL, NULL);

      if (!input_x || !output)
      {
          printf("Error: Failed to allocate device memory!\n");
          exit(1);
      }

      // Write our data set into the input array in device memory 
      //
      auto start = high_resolution_clock::now();
      err = clEnqueueWriteBuffer(commands, input_x, CL_TRUE, 0, sizeof(G1<mnt4753_pp>), data_x, 0, NULL, NULL);
      if (err != CL_SUCCESS)
      {
          printf("Error: Failed to write to source array!\n");
          exit(1);
      }
      err = clEnqueueWriteBuffer(commands, input_y, CL_TRUE, 0, sizeof(G1<mnt4753_pp>) * count, data_y, 0, NULL, NULL);
      if (err != CL_SUCCESS)
      {
          printf("Error: Failed to write to source array!\n");
          exit(1);
      }
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<microseconds>(stop - start); 
      cout << "Time taken by GPU write function: "
        << duration.count() << " microseconds" << endl;

      // Set the arguments to our compute kernel
      //
      err = 0;
      err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &input_x);
      err  = clSetKernelArg(kernel, 1, sizeof(cl_mem), &input_y);
      err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &output);
      err |= clSetKernelArg(kernel, 3, sizeof(unsigned int), &count);
      if (err != CL_SUCCESS)
      {
          printf("Error: Failed to set kernel arguments! %d\n", err);
          exit(1);
      }

      // Get the maximum work group size for executing the kernel on the device
      //
      err = clGetKernelWorkGroupInfo(kernel, devices[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(local), &local, NULL);
      if (err != CL_SUCCESS)
      {
          printf("Error: Failed to retrieve kernel work group info! %d\n", err);
          exit(1);
      }

      printf("Max work size: %u\n", local);

      // Execute the kernel over the entire range of our 1d input data set
      // using the maximum number of work group items for this device
      //
      global = count;
      //global = 1;
      printf("queueing kernel\n");
      err = clEnqueueNDRangeKernel(commands, kernel, 1, NULL, &global, &local, 0, NULL, &event);
      if (err)
      {
          printf("Error: Failed to execute kernel!\n");
          return EXIT_FAILURE;
      }

      clWaitForEvents(1, &event);
      clFinish(commands);

      // Time kernel execution time without read/write
      //
      clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
      clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

      double nanoSeconds = time_end-time_start;
      printf("OpenCl Execution time is: %0.3f milliseconds \n",nanoSeconds / 1000000.0);

      // Read back the results from the device to verify the output
      //
      start = high_resolution_clock::now();
      err = clEnqueueReadBuffer( commands, output, CL_TRUE, 0, sizeof(G1<mnt4753_pp>), results, 0, NULL, NULL );  
      if (err != CL_SUCCESS)
      {
          printf("Error: Failed to read output array! %d\n", err);
          exit(1);
      }
      stop = high_resolution_clock::now();
      duration = duration_cast<microseconds>(stop - start); 
      cout << "Time taken by GPU read function: "
        << duration.count() << " microseconds" << endl;
      // Validate our results
      //
      printf("Kernel Result \n");
      results[0].print();

      // results[0].coeff_a.mont_repr.print_hex();
      // for(int i=0; i<12; i++) {
      //   //std::cout << "Length of array = " << (sizeof(results[1013].non_residue.mont_repr.data)/sizeof(*results[1013].non_residue.mont_repr.data)) << std::endl;
      //   cl_uint x;
      //   cl_uint y;
      //   x = (cl_uint)((results[0].coeff_a.mont_repr.data[i] & 0xFFFFFFFF00000000LL) >> 32);
      //   y = (cl_uint)(results[0].coeff_a.mont_repr.data[i] & 0xFFFFFFFFLL);
      //   gmp_printf("%Mx\n", results[0].coeff_a.mont_repr.data[i]);
      //   printf("%x\n", x);
      //   printf("%x\n", y);
      // }

      results[0].zero().print_coordinates();

      // for(int i=0; i<12; i++) {
      //   //printf("%x\n", results[1013].c0.mod.data[i]);
      //   //std::cout << "Length of array = " << (sizeof(results[1013].non_residue.mont_repr.data)/sizeof(*results[1013].non_residue.mont_repr.data)) << std::endl;
      //   cl_uint x;
      //   cl_uint y;
      //   x = (cl_uint)((results[1013].c0.one().mont_repr.data[i] & 0xFFFFFFFF00000000LL) >> 32);
      //   y = (cl_uint)(results[1013].c0.one().mont_repr.data[i] & 0xFFFFFFFFLL);
      //   gmp_printf("%Mx\n", results[1013].c0.one().mont_repr.data[i]);
      //   printf("%x\n", x);
      //   printf("%x\n", y);
      // }

      printf("CPU Result\n");
      G1<mnt4753_pp> _h4_1 = G1<mnt4753_pp>::zero();

      //for (size_t i = 0; i < n; ++i) { _h4_1 = _h4_1 + g4_1[i]; }
      _h4_1 = _h4_1.add(g4_1[0]);
      _h4_1 = _h4_1.add(g4_1[1]);
      _h4_1 = g4_1[6].add(g4_1[6]);
      _h4_1 = _h4_1.add(g4_1[3]);
      _h4_1 = g4_1[4].add(g4_1[4]);
      _h4_1.print();
      //  g4_1[1].X().print();
      correct = 0;

      // there is some fuckery on the results fqe struct, cant equality check mont_repr
      if(results[0] == _h4_1) {
        correct++;
      }

      
      // Print a brief summary detailing the results
      //
      printf("Computed '%d/%d' correct fq3 values!\n", correct, count);
      //x0[1014].mont_repr.print();
      //x1[1014].mont_repr.print();
      // Shutdown and cleanup
      //
      clReleaseMemObject(input_x);
      clReleaseMemObject(input_y);
      clReleaseMemObject(output);
      clReleaseProgram(program);
      clReleaseKernel(kernel);
      clReleaseCommandQueue(commands);
      clReleaseContext(context);

      // OPENCL END
      //break;
    }
    fclose(outputs);

    return 0;
}
